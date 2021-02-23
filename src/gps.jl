# Methods related to all things GP regression

struct GPInfo
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
Generate GP regression given the dataset.
"""
function generate_estimates(params, x_train, y_train; reuse_gp_flag=false, filename_appendix=nothing)
    exp_dir = create_experiment_directory(params)
    # Generate a set of GPRs if none are provided.
    gps_dir = @sprintf("%s/gps", exp_dir)
    !isdir(gps_dir) && mkpath(gps_dir)
    gps_filename = @sprintf("%s-m%d-σ%1.3f-rs%d-gps", params.system_params.mode_tag, params.data_params.data_num, params.data_params.noise_sigma, params.random_seed)
    gps_filename = isnothing(filename_appendix) ? gps_filename : "$gps_filename-$filename_appendix"
    local_gp_file = "$exp_dir/$gps_filename.bin"
    data_deps = params.system_params.dependency_dims

    if isfile(local_gp_file) && reuse_gp_flag
        @info "Resuing gp: " local_gp_file 
        open(local_gp_file) do f
            gp_set = deserialize(f)
        end
    else 
        # Train GPs on this system
        gp_set = Dict()
        # TODO: This is not hyperparameterized.
        ls = 0.65
        for (i, out_dim) in enumerate(keys(data_deps)) 
            # Handle data dependency here
            x_train_sub = x_train[:, findall(.>(0), data_deps[out_dim])[:]]
            m_prior = MeanZero()
            k_prior = SE(ls, 0.)
            lnoise = log(sqrt(1+2/length(x_train_sub))) # Generalize to handle any bound
            # opt_idx = StatsBase.sample(1:length(y_train[:,1]), params.data_params.data_num_optimize, replace = false)
            # gp_pre = GP(x_train_sub[opt_idx, :]', y_train[opt_idx,i], m_prior, k_prior, lnoise) 
            # optimize!(gp_pre)
            gp = GP(x_train_sub', y_train[:,i], m_prior, k_prior, lnoise)
            gp_set["x$i"] = deepcopy(gp)
        end
    end

    gp_info_dict = Dict()
    for dim_key in keys(data_deps) 
        gp_info_dict[dim_key] = create_gp_info(params, gp_set, dim_key) 
    end

    # TODO: Plotting should be exposed to get accurate timing
    # create_plots ? plot_gp_fields(exp_dir, dyn_fn) : nothing

    return gp_set, gp_info_dict
end

"""
Save the GP info and objects.
"""
function save_gp_info(params, gp_set, gp_info_dict; save_global_gps=true, filename_appendix=nothing)
    exp_dir = create_experiment_directory(params)
    if save_global_gps
        gps_dir = @sprintf("%s/gps", params.experiment_directory)
        !isdir(gps_dir) && mkpath(gps_dir)
        gps_filename = @sprintf("%s/%s-m%d-σ%1.3f-rs%d-gps.bin", gps_dir, params.system_params.mode_tag, params.data_params.data_num, params.data_params.noise_sigma, params.random_seed)
        # Save the GPs for further analysis. 
        @info "Saving GP regressions to experiment directory..."
        open(gps_filename, "w") do f
            serialize(f, gp_set)
        end
    end
   
    gp_save_info = Dict()
    gp_save_info[:gp_set] = gp_set
    gp_save_info[:gp_info] = gp_info_dict

    local_gp_dir = "$exp_dir/gps"
    !isdir(local_gp_dir) && mkpath(local_gp_dir)
    gps_filename = @sprintf("%s-m%d-σ%1.3f-rs%d-gps", params.system_params.mode_tag, params.data_params.data_num, params.data_params.noise_sigma, params.random_seed)
    gps_filename = isnothing(filename_appendix) ? gps_filename : "$gps_filename-$filename_appendix"
    local_gp_file = "$local_gp_dir/$gps_filename.bin"
    open(local_gp_file, "w") do f
        serialize(f, gp_save_info)
    end
end

"""
Update the GP with the datapoint given as an input-output tuple. 
"""
function update_gp(gp, datapoint::Tuple)
    # Update the GP in place - does not return a new GP structure. 
    newx = [gp.x'; datapoint[1][:]']'
    newy = gp.y
    push!(newy, datapoint[2])
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
