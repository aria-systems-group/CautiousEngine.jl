" Generate image overapproximations for each discrete state. 
# Arguments
- `params::ExperimentParameters` - Experiment parameters structure
- `gp_info_dict::Dict` - Dictionary with GP set and RKHS info
"
function generate_region_images(params::ExperimentParameters, gp_info_dict::Dict)
	region_dict = create_region_data(params.domain, params.discretization_step)
	region_post_dict = Dict()

	# Minimum and maximum extents
	domain = params.domain
    # TODO: put this in a different fcn?
	min_ex = [domain[dim_key][1] for dim_key in keys(domain)]
	max_ex = [domain[dim_key][2] for dim_key in keys(domain)]
	mx = zeros(2,1)
	mx[1] = max_ex[1]*1.1
	mx[2] = max_ex[2]*1.1
	# Store the predicted mean and covariance for each sampled point 
	dim_keys = keys(domain)
	σ_ubs = Dict()

	# Calculate GP UBs
	for dim_key in dim_keys
		σ_ubs[dim_key] = 1.1*sqrt(predict_f(gp_info_dict[dim_key].gp, mx)[2][1])
	end

    # TODO: Handle this flag better
	# σ_ubs = nothing
	
    Threads.@threads for i=1:region_dict.count-1
        region_post_dict[i] = bound_extent(region_dict[i], gp_info_dict, params.system_params.dependency_dims; known_component=params.system_params.known_dynamics_fcn, σ_ubs=σ_ubs)
    end

	region_data = Dict(:extents=>region_dict, :posts=>region_post_dict)

    return region_data
end

" Generate image overapproximations for each discrete state using local regressions. 
# Arguments
- `params::ExperimentParameters` - Experiment parameters structure
- `x_train::AbstractArray{Any}` - Array of input data 
- `y_train::AbstractArray{Any}` - Array of input data 
"
function generate_region_images(params::ExperimentParameters, x_train::AbstractArray{Any}, y_train::AbstractArray{Any}; reuse_regions_flag=false)

    exp_dir = params.experiment_directory
    region_filename = "$exp_dir/regions.bson"

    # TODO: Remove this - reloading is handled elsewhere
    if reuse_regions_flag && isfile(region_filename)
        # TODO: Replace this fella!!!!
        region_data = BSON.load(region_filename)
    else
        region_dict = create_region_data(params.domain, params.discretization_step)
        num_regions = region_dict.count 
        region_post_dict = Dict()
        region_gp_dict = Dict()

        # Minimum and maximum extents
        domain = params.domain
        min_ex = [domain[dim_key][1] for dim_key in keys(domain)]
        max_ex = [domain[dim_key][2] for dim_key in keys(domain)]

        # Store the predicted mean and covariance for each sampled point 
        dim_keys = keys(domain)
        Threads.@threads for i=1:num_regions
            # Get subset of data here
            if i == num_regions
                extent = region_dict[i]
                lb = [extent[dim_key][1] for dim_key in dim_keys]
                ub = [extent[dim_key][2] for dim_key in dim_keys]
                # Assume 2D
                gp_set, gp_info_dict = generate_estimates(params, x_train, y_train, filename_appendix=i) 
                save_gp_info(params, gp_set, gp_info_dict, save_global_gps=false, filename_appendix=i)
                region_post_dict[i] = bound_extent(extent, lb, ub, gp_info_dict, dim_keys, params.system_params.dependency_dims; known_part_flag=!isnothing(params.system_params.known_dynamics_fcn))
                region_gp_dict[i] = gp_info_dict 
            else
                extent = region_dict[i]
                lb = [extent[dim_key][1] for dim_key in dim_keys]
                ub = [extent[dim_key][2] for dim_key in dim_keys]
                # Assume 2D
                center = [mean(lb), mean(ub)]
               
                kdtree = KDTree(x_train);
                
                num_neighbors = params.data_params.num_neighbors
                sub_idx, _ = knn(kdtree, center, num_neighbors, true)
                
                x_sub = @view x_train[:, sub_idx]
                y_sub = @view y_train[:, sub_idx]

                gp_set, gp_info_dict = generate_estimates(params, x_sub, y_sub, filename_appendix=i) 
                # TODO: saving GPs might be a dumb idea to get a true sense of time.
                save_gp_info(params, gp_set, gp_info_dict, save_global_gps=false, filename_appendix=i)
                region_post_dict[i] = bound_extent(extent, lb, ub, gp_info_dict, dim_keys, params.system_params.dependency_dims; known_part_flag=!isnothing(params.system_params.known_dynamics_fcn))
                region_gp_dict[i] = gp_info_dict
            end
        end  # End threaded forloop

        region_data = Dict()
        region_data[:extents] = region_dict
        region_data[:posts] = region_post_dict
        region_data[:gps] = region_gp_dict
    end

    return region_data
end

" Generate overapproximations of posterior mean and covariance functions using one of several methods.
# Arguments
- `extent::Dict` - Discrete state extent 
- `gp_info_dict::Dict` - Dictionary with GP set and RKHS info
- `data_deps::Dict` - Dictionary indicating input-output data dependencies
"
function bound_extent(extent::Dict, gp_info_dict::Dict, data_deps; known_component=nothing, σ_ubs=nothing, σ_approx_flag=false)
    # TODO: mod keyword handling and document
    dim_keys = keys(extent)
    lb = [extent[dim_key][1] for dim_key in dim_keys]
    ub = [extent[dim_key][2] for dim_key in dim_keys]
    post_extent = Dict()
    σ_bounds = []

    for dim_key in dim_keys 
        lbf = lb[findall(.>(0), data_deps[dim_key][:])]
        ubf = ub[findall(.>(0), data_deps[dim_key][:])]
        x_lb, μ_L_lb, μ_L_ub = compute_μ_bounds_bnb(deepcopy(gp_info_dict[dim_key].gp), lbf, ubf) 
        x_ub, μ_U_lb, μ_U_ub = compute_μ_bounds_bnb(deepcopy(gp_info_dict[dim_key].gp), lbf, ubf, max_flag=true)

        if σ_approx_flag
            _, σ_U_lb, σ_U_ub = compute_σ_ub_bounds_approx(gp_info_dict[dim_key].gp, lbf, ubf) 
        elseif !isnothing(σ_ubs)
            _, σ_U_lb, σ_U_ub = compute_σ_ub_bounds_from_gp(gp_info_dict[dim_key].gp, lbf, ubf, ub=σ_ubs[dim_key])
        else
            _, σ_U_lb, σ_U_ub = compute_σ_ub_bounds(gp_info_dict[dim_key].gp, gp_info_dict[dim_key].Kinv, lbf, ubf)
        end
       
        if !(μ_L_lb <= μ_U_ub)
            @error "μ LB is greater than μ UB! ", μ_L_lb, μ_U_ub
            throw("aborting") 
            @assert μ_L_lb <= μ_U_ub
        end
        if !(σ_U_lb <= σ_U_ub)
            @error "Sigma LB is greater than sigma UB! ", σ_U_lb, σ_U_ub
            throw("aborting")
            @assert σ_U_lb <= σ_U_ub
        end 
        if !isnothing(known_component) 
            # TODO: Assume an identity form for now. Requires its own bounding process.
            post_extent[dim_key] = [extent[dim_key][1] + μ_L_lb, extent[dim_key][2] + μ_U_ub] 
        else
            post_extent[dim_key] = [μ_L_lb, μ_U_ub]
        end
        push!(σ_bounds, σ_U_ub)
    end
    return post_extent, σ_bounds
end