# Methods related to all things GP regression

mutable struct GPInfo
    gp
    γ_bound::Float64
    RKHS_bound::Float64
    bound_type::String
    post_scale_factor::Float64
    Kinv
end

"""
Select datapoints in a neighborhood of a point.
"""
function get_points_in_neighborhood(center, radius, x_train, y_train)
    norm_array = [norm(row) for row in eachrow(center .- x_train)]
    bool_array = norm_array .<= radius
    x_sub = @view x_train[bool_array, :]
    y_sub = @view y_train[bool_array, :]
    return x_sub, y_sub
end

"""
Create GP From Given Data
"""
function train_gps(x_train, y_train; se_params=[0., 0.65], optimize_hyperparameters=false, lnoise=nothing, opt_fraction=1.0)
    gp_set = Dict()
    dim_keys = ["x$i" for i=1:size(x_train,1)]
    for (i, out_dim) in enumerate(dim_keys) 
        # Handle data dependency here
        # x_train_sub = x_train[:, findall(.>(0), data_deps[out_dim])[:]]
        x_train_sub = x_train
        gp = train_gp_1dim(x_train, y_train[i,:]; se_params=se_params, optimize_hyperparameters=optimize_hyperparameters, lnoise=lnoise, opt_fraction=opt_fraction)
        
        gp_set["x$i"] = deepcopy(gp)
    end
    return gp_set
end

function generate_gp_filename(params; filename_appendix=nothing)
    gps_dir = @sprintf("%s/gps", params.experiment_directory)
    !isdir(gps_dir) && mkpath(gps_dir)
    gps_filename = @sprintf("%s-m%d-σ%1.3f-rs%d-gps", params.system_params.mode_tag, params.data_params.data_num, params.data_params.noise_sigma, params.random_seed)
    gps_filename = isnothing(filename_appendix) ? gps_filename : "$gps_filename-$filename_appendix"
    return "$gps_dir/$gps_filename.bin" 
end

function train_gp_1dim(x_train, y_train; se_params=[0., 0.65], optimize_hyperparameters=false, lnoise=nothing, opt_fraction=1.0)
    m_prior = MeanZero()
    k_prior = SE(se_params[2], se_params[1])

    if isnothing(lnoise)
        lnoise = log(sqrt(1+2/size(x_train, 2))) # Generalize to handle any bound
    end

    if optimize_hyperparameters
        @assert 0 < opt_fraction <= 1
        num_opt = Int(opt_fraction*length(x_train))
        opt_idx = StatsBase.sample(1:length(y_train), num_opt, replace = false)
        gp_pre = GP(x_train[:, opt_idx], y_train[opt_idx], m_prior, k_prior, lnoise) 
        optimize!(gp_pre)
        gp = update_gp_data(gp_pre, x_train, y_train)
    else
        gp = GP(x_train, y_train, m_prior, k_prior, lnoise)
        if optimize_hyperparameters
            optimize!(gp)
        end 
    end

    return gp
end
"""
Generate GP regression given the dataset.
"""
function generate_estimates(params, x_train, y_train)
    data_deps = params.system_params.dependency_dims
    x_train_sub = x_train #[:, findall(.>(0), data_deps[out_dim])[:]]
    optimize_hyperparameters = params.data_params.data_frac_optimize > 0. ? true : false
    gp_set = train_gps(x_train, y_train; se_params=[0., 0.65], optimize_hyperparameters=optimize_hyperparameters, lnoise=nothing, opt_fraction=params.data_params.data_frac_optimize)

    gp_info_dict = Dict()
    for dim_key in keys(data_deps) 
        gp_info_dict[dim_key] = create_gp_info(params, gp_set, dim_key) 
    end

    return gp_set, gp_info_dict
end

"""
Create Local GP From Global GP
"""
function create_local_gp_info_from_global(x_data, y_data, global_gp_info_dict)

    local_gp_info_dict = copy(global_gp_info_dict) 

    for (i, gp_key) in enumerate(keys(global_gp_info_dict))
        gp = global_gp_info_dict[gp_key].gp
        gp_local = update_gp_data(gp, x_data, y_data)
        local_gp_info_dict[gp_key].gp = gp_local
        update_gp_info(local_gp_info_dict[gp_key])
    end

    return local_gp_info_dict
end

"""
Create the GP Info structure componenets.
"""
function create_gp_info(params, gp_set, dim_key)
    gp = gp_set[dim_key]

    scale_factor = params.data_params.bound_type == "rkhs-tight" ? params.data_params.noise_sigma/sqrt(1. + 2. / gp.nobs) : 1.
    domain = params.domain
    diam_domain = 0.
    for dim_key in keys(domain)
        diam_domain += (domain[dim_key][1] - domain[dim_key][2])^2
    end
    diam_domain = sqrt(diam_domain)
    σ_inf = sqrt(gp.kernel.σ2*exp(-1/2*(diam_domain)^2/gp.kernel.ℓ2))

    # Calculating the RKHS parameter bounds
    RKHS_bound = abs(domain[dim_key][2] + params.system_params.lipschitz_bound*diam_domain)/ σ_inf
    B = 1 + (1 + 2/(gp.nobs))^(-1)
    γ = 0.5*gp.nobs*log(B)
    K_inv = inv(gp.cK.mat + exp(gp.logNoise.value)^2*I)
    gp_info = GPInfo(gp, γ, RKHS_bound, params.data_params.bound_type, scale_factor, K_inv)

    return gp_info
end

"""
Predict mean function
"""
function predict_N_steps(gp_info_dict, x0, N)
    xp_trajectory = []
    xc = x0
    dim_keys = keys(gp_info_dict)

    for i=1:N
        xp = zeros(length(dim_keys), 1) 
        for (j, key) in enumerate(dim_keys)
            μ, _ = predict_f(gp_info_dict[key].gp, xc)
            xp[j] = μ[1] 
        end

        push!(xp_trajectory, xp)
        xc = xp
    end

    return xp_trajectory
end

"""
Predict mean and vairance over one step
"""
function predict_one_step_full(gp_info_dict, x0; known_flag=true)
    μ_dict = Dict()
    σ_dict = Dict()
    dim_keys = keys(gp_info_dict)
    for (i, key) in enumerate(dim_keys)
        μ, σ = predict_f(gp_info_dict[key].gp, x0)
        μ_dict[key] = μ[1] + x0[i]
        σ_dict[key] = σ[1]
    end
return μ_dict, σ_dict
end

"""
Update the GP Info after updating the gp
"""
function update_gp_info(gp_info)
    gp = gp_info.gp 
    bound_type = gp_info.bound_type

    old_nobs = size(gp_info.Kinv)[1]
    scale_factor = gp_info.bound_type == "rkhs-tight" ? gp_info.post_scale_factor*sqrt(1. + 2. / (old_nobs))/sqrt(1. + 2. / gp.nobs) : 1.

    B = 1 + (1 + 2/(gp.nobs))^(-1)
    γ = 0.5*gp.nobs*log(B)

    K_inv = inv(gp.cK.mat + exp(gp.logNoise.value)^2*I)

    gp_info.post_scale_factor = scale_factor
    gp_info.Kinv = K_inv
    gp_info.γ_bound = γ

    return nothing
end

"""
Save the GP info and objects.
"""
function save_gp_info(params, gp_set, gp_info_dict; filename_appendix=nothing)
    exp_dir = params.experiment_directory 

    gp_save_info = Dict()
    gp_save_info[:gp_set] = gp_set
    gp_save_info[:gp_info] = gp_info_dict

    gps_filename = generate_gp_filename(params, filename_appendix=filename_appendix)

    # Save the GPs for further analysis. 
    @info "Saving GP regressions to experiment directory..."
    open(gps_filename, "w") do f
        serialize(f, gp_save_info)
    end
end

function load_gp_info(filename)
    gp_save_info = nothing
    open(filename, "r") do f
        gp_save_info = deserialize(f)
    end
    return gp_save_info[:gp_set], gp_save_info[:gp_info]
end

"""
Update the GP with the specified dataset and return a new GP object
"""
function update_gp_data(gp, x_data, y_data)
    new_gp = deepcopy(gp)
    GaussianProcesses.fit!(new_gp, x_data, y_data)
    return new_gp
end


"""
Update the GP by adding a single datapoint
"""
function add_datapoint_to_gp(gp::GPE, datapoint::Tuple)
    newx = [gp.x datapoint[1]]
    newy = [gp.y; datapoint[2]]
    GaussianProcesses.fit!(gp, newx, newy)
end

"""
Load a GP bin file.
"""
function load_gp(filename)
    f = open(filename)
    gp_set = deserialize(f)
    close(f)
    return gp_set
end
