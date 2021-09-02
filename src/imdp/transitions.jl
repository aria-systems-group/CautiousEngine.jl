
"""
Generates transition bounds with the given parameters.
"""
function generate_transition_bounds(params, gp_info_dict, region_data; 
                                    reuse_mats_flag=false) 

    exp_dir = params.experiment_directory
    region_dict = region_data[:extents]
    region_post_dict = region_data[:posts]

    tmat_filename = "$exp_dir/transition_mats.bson"
    if !reuse_mats_flag || !isfile(tmat_filename)
        @info "Calculating transition probability bounds between regions..."
        minPr_mat, maxPr_mat = process(region_dict, region_post_dict, params.data_params.epsilon, gp_info_dict, params)
    end

    res_mats = Dict("minPr" => minPr_mat, "maxPr" => maxPr_mat)
    return res_mats 
end

"""
Generates transition bounds with the given parameters.
"""
function generate_transition_bounds(params, region_data; 
                                    reuse_mats_flag=false) 

    exp_dir = params.experiment_directory
    region_dict = region_data[:extents]
    region_post_dict = region_data[:posts]
    region_gp_dict = region_data[:gps]

    tmat_filename = "$exp_dir/transition_mats.bson"
    if !reuse_mats_flag || !isfile(tmat_filename)
        @info "Calculating transition probability bounds between regions..."
        results_df = process_foo(region_dict, region_post_dict, region_gp_dict, params.data_params.epsilon, params)

        # Save the matrices for use in MATLAB verification tool
        num_states = region_dict.count
        minPr_mat = spzeros(num_states, num_states)
        maxPr_mat = spzeros(num_states, num_states)
        for i = 1:1:num_states
            sub_row = results_df[results_df.Set1 .== i, :]
            if i == num_states
                minPr_mat[i, i] = 1.
                maxPr_mat[i, i] = 1.
            else
                for j = 1:1:num_states
                    if j == num_states
                        subsub_row = sub_row[sub_row.Set2 .== j, :]
                        minPr_mat[i, j] = 1. - subsub_row.MaxPr[1]
                        maxPr_mat[i, j] = 1. - subsub_row.MinPr[1]
                    else
                        subsub_row = sub_row[sub_row.Set2 .== j, :]
                        minPr_mat[i, j] = subsub_row.MinPr[1]
                        maxPr_mat[i, j] = subsub_row.MaxPr[1]
                    end
                end
            end
        end
    end

    res_mats = Dict("minPr" => minPr_mat, "maxPr" => maxPr_mat)
    return res_mats 
end

function load_transition_matrices(filename)
    res_mats = BSON.load(filename)
    return res_mats
end

function save_transition_matrices(params::ExperimentParameters, res_mats::Dict)
    exp_dir = params.experiment_directory
    tmat_filename = "$exp_dir/transition_mats.bson"
    save_transition_matrices(res_mats, tmat_filename) 
end

function save_transition_matrices(res_mats::Dict, filename::String)
    bson(filename, res_mats) 
end

function fetch_results(N, r)
    results_df = DataFrame(Set1 = Int[], Set2 = Int[], MinPr = Float64[], MaxPr = Float64[], MinPrPoint = Array[], MaxPrPoint = Array[])
    for j in 1:N
        total_results = fetch(r[j])
        push!(results_df, total_results)
    end
    return results_df 
end

function process(region_dict, region_post_dict, epsilon, gp_dict, params)

    minPrs = spzeros(region_dict.count, region_dict.count)
    maxPrs = spzeros(region_dict.count, region_dict.count)

    mean_dict = Dict()
    for i=1:region_dict.count-1
        region = region_dict[i]
        mean_dict[i] = [mean(region[dim_key]) for dim_key in keys(region)] 
    end

    set_radius = sqrt(2)*params.discretization_step["x1"]

    for i=1:region_dict.count-1
        ϵ_crit = calculate_ϵ_crit(gp_dict, region_post_dict[i][2])
        # calculate center of post image
        post_image = region_post_dict[i][1]
        mean_pt = [mean(post_image[dim_key]) for dim_key in keys(post_image)]
        # TODO: This is not generalized
        image_radius = sqrt((mean_pt[1] - post_image["x1"][1])^2 + (mean_pt[1] - post_image["x2"][1])^2) 
        for j=1:region_dict.count-1
            if fast_check(mean_pt, mean_dict[j], ϵ_crit, image_radius, set_radius)
                res = region_pair_transitions((i,j), region_dict, region_post_dict, epsilon, gp_dict, params)
            else
                res = [i, j, 0., 0., [-1., -1.], [-11., -11.]]
            end
            minPrs[i, j] = res[3]
            maxPrs[i, j] = res[4]
        end

        res = region_pair_transitions((i,region_dict.count), region_dict, region_post_dict, epsilon, gp_dict, params)
        minPrs[i, region_dict.count] = 1. - res[4]
        maxPrs[i, region_dict.count] = 1. - res[3]
    end

    minPrs[region_dict.count, :] .= 1.
    maxPrs[region_dict.count, :] .= 1.
   
    return minPrs, maxPrs 
end

function process_foo(region_dict, region_post_dict, region_gp_dict, epsilon, params)
    N = region_dict.count
    r = Array{Array, 1}(undef, N)

    for (i, pair) = enumerate(Base.product(1:N, 1:N))
        i == N ? continue : nothing
        r[i] = region_pair_transitions(pairs[i], region_dict, region_post_dict, epsilon, region_gp_dict[pairs[i][1]], params)
    end

    results_df = fetch_results(N, r)
    return results_df
end

"Determine the set of states that will have non-zero transition probability upper bounds."
function fast_check(mean_pt, mean_target, ϵ_crit, image_radius, set_radius)
    flag = norm(mean_pt - mean_target) < ϵ_crit + image_radius + set_radius
    return flag
end

"Calculate the transition probabilities between a pair of regions."
function region_pair_transitions(region_pair, region_dict, region_post_dict, ϵ, gp_info_dict, params)
    dict_row = region_post_dict[region_pair[1]]
    σ_bounds = dict_row[2] 

    # Get the subset of the safe region that matters.
    safe_keys = params.data_params.safety_dims
    S = Dict()
    Q_post_extent = Dict()
    for dim_key in safe_keys
        S[dim_key] = region_dict[region_pair[2]][dim_key] 
        Q_post_extent[dim_key] = dict_row[1][dim_key]
    end

    # For now, do a naiive shrink - first noise, then the next
    if !isnothing(params.system_params.process_noise_dist)
        ϵ_noise = params.data_params.eta
        Pr_noise = (cdf(params.system_params.process_noise_dist, ϵ_noise) - cdf(params.system_params.process_noise_dist, -ϵ_noise))^length(σ_bounds)
    else
        ϵ_noise = 0.
        Pr_noise = 1.
    end

    Sshrink = margin_extent(S, ϵ_noise)
    Sexpand = margin_extent(S, -ϵ_noise)
    s_shape = extent_to_shape(Sshrink)
    S_shape = extent_to_shape(Sexpand)

    Q_shape = extent_to_shape(Q_post_extent)
    n_dims = length(safe_keys)
    axes = I+zeros(n_dims, n_dims)

    if ϵ > 0
        # Manual Epsilon Choice
        sd_min = Inf
        for dim_key in keys(Sshrink)
            sd_min = minimum([sd_min, 0.5*abs(Sshrink[dim_key][1]-Sshrink[dim_key][2])])
        end
        shrink_ϵ = minimum([ϵ, sd_min])
        expand_ϵ = ϵ

        Sshrink2 = margin_extent(S, shrink_ϵ)
        Sexpand2 = margin_extent(S, -expand_ϵ)
        s_shape2 = extent_to_shape(Sshrink2)
        S_shape2 = extent_to_shape(Sexpand2)

        shrink_int, _ = does_shape2_contain_shape1(Q_shape, s_shape2, axes=axes)
        expand_int, _ = do_shapes_intersect(Q_shape, S_shape2, axes=axes)
    else
        # Auto Epsilon Choice
        shrink_int, shrink_ϵ = does_shape2_contain_shape1(Q_shape, s_shape, axes=axes)
        expand_int, expand_ϵ = do_shapes_intersect(Q_shape, S_shape, axes=axes)
    end

    if (shrink_int && !expand_int)
        @info "Shrink Flag: ", shrink_int, "Expand Flag: ", expand_int
        @info "Q Post Extent: ", Q_post_extent 
        @info "Q Post Shape: ", Q_shape
        @info "S Extent: ", S
        @info "S Shrink Ex: ", Sshrink
        @info "S Shrink Shape: ", s_shape
        @info "S Expand Ex: ", Sexpand
        @info "S Expand Shape: " S_shape
        @info "Shrink Ep: ", shrink_ϵ, "Expand Ep: ", expand_ϵ 
        @error "You cannot both interesect and not intersect at the same time!"
    end

    if shrink_int 
        P_bound_lower = calculate_probability_bound(gp_info_dict, σ_bounds, shrink_ϵ)
        minPr_LB = P_bound_lower*Pr_noise 
    else
        minPr_LB = 0.
    end

    if !expand_int 
        Pr_bound_upper = calculate_probability_bound(gp_info_dict, σ_bounds, expand_ϵ)
        maxPr_UB = 1. - Pr_bound_upper
    else
        maxPr_UB = 1.
    end

    @assert maxPr_UB >= minPr_LB

    foo = [0., 0.]
    frame_row = [region_pair[1], region_pair[2], minPr_LB, maxPr_UB, [shrink_ϵ, expand_ϵ], foo]
    return frame_row
end

""" 
Calculate transition probability bounds between two regions with a different flavor
"""
function region_pair_transitions(region_pair, region_dict, region_post_dict, ϵ, gp_info_dict, process_noise_dist, ϵ_noise)
    dict_row = region_post_dict[region_pair[1]]
    σ_bounds = dict_row[2] 

    # Get the subset of the safe region that matters.
    # safe_keys = params.data_params.safety_dims
    S = Dict()
    Q_post_extent = Dict()
    safe_keys = keys(dict_row[1])
    for dim_key in safe_keys
        S[dim_key] = region_dict[region_pair[2]][dim_key] 
        Q_post_extent[dim_key] = dict_row[1][dim_key]
    end

    # For now, do a naiive shrink - first noise, then the next
    if !isnothing(process_noise_dist)
        Pr_noise = (cdf(process_noise_dist, ϵ_noise) - cdf(process_noise_dist, -ϵ_noise))^length(σ_bounds)
    else
        ϵ_noise = 0.
        Pr_noise = 1.
    end

    Sshrink = margin_extent(S, ϵ_noise)
    Sexpand = margin_extent(S, -ϵ_noise)
    s_shape = extent_to_shape(Sshrink)
    S_shape = extent_to_shape(Sexpand)

    Q_shape = extent_to_shape(Q_post_extent)
    n_dims = length(safe_keys)
    axes = I+zeros(n_dims, n_dims)

    if ϵ > 0
        # Manual Epsilon Choice
        sd_min = Inf
        for dim_key in keys(Sshrink)
            sd_min = minimum([sd_min, 0.5*abs(Sshrink[dim_key][1]-Sshrink[dim_key][2])])
        end
        shrink_ϵ = minimum([ϵ, sd_min])
        expand_ϵ = ϵ

        Sshrink2 = margin_extent(S, shrink_ϵ)
        Sexpand2 = margin_extent(S, -expand_ϵ)
        s_shape2 = extent_to_shape(Sshrink2)
        S_shape2 = extent_to_shape(Sexpand2)

        shrink_int, _ = does_shape2_contain_shape1(Q_shape, s_shape2, axes=axes)
        expand_int, _ = do_shapes_intersect(Q_shape, S_shape2, axes=axes)
    else
        # Auto Epsilon Choice
        shrink_int, shrink_ϵ = does_shape2_contain_shape1(Q_shape, s_shape, axes=axes)
        expand_int, expand_ϵ = do_shapes_intersect(Q_shape, S_shape, axes=axes)
    end

    if (shrink_int && !expand_int)
        @info "Shrink Flag: ", shrink_int, "Expand Flag: ", expand_int
        @info "Q Post Extent: ", Q_post_extent 
        @info "Q Post Shape: ", Q_shape
        @info "S Extent: ", S
        @info "S Shrink Ex: ", Sshrink
        @info "S Shrink Shape: ", s_shape
        @info "S Expand Ex: ", Sexpand
        @info "S Expand Shape: " S_shape
        @info "Shrink Ep: ", shrink_ϵ, "Expand Ep: ", expand_ϵ 
        @error "You cannot both interesect and not intersect at the same time!"
    end

    if shrink_int 
        P_bound_lower = calculate_probability_bound(gp_info_dict, σ_bounds, shrink_ϵ)
        minPr_LB = P_bound_lower*Pr_noise 
    else
        minPr_LB = 0.
    end

    if !expand_int 
        Pr_bound_upper = calculate_probability_bound(gp_info_dict, σ_bounds, expand_ϵ)
        maxPr_UB = 1. - Pr_bound_upper
    else
        maxPr_UB = 1.
    end

    @assert maxPr_UB >= minPr_LB

    foo = [0., 0.]
    frame_row = [region_pair[1], region_pair[2], minPr_LB, maxPr_UB, [shrink_ϵ, expand_ϵ], foo]
    return frame_row
end

"""
Calculate the transition probability bounds from a point to a target region.
"""
function point_region_transitions(μ::Dict, σ::Dict, target_region_idx, region_dict, ϵ, gp_info_dict; process_noise_dist=nothing)

    # Get the subset of the safe region that matters.
    # safe_keys = params.data_params.safety_dims
    safe_keys = keys(gp_info_dict)
    S = Dict()
    Q_post_extent = Dict()
    for dim_key in safe_keys
        S[dim_key] = region_dict[target_region_idx][dim_key] 
        Q_post_extent[dim_key] = [μ[dim_key], μ[dim_key]] 
    end

    n_dims = length(μ)

    # For now, do a naiive shrink - first noise, then the next
    if !isnothing(process_noise_dist)
        ϵ_noise = 0.02
        Pr_noise = (cdf(process_noise_dist, ϵ_noise) - cdf(process_noise_dist, -ϵ_noise))^n_dims
    else
        ϵ_noise = 0.
        Pr_noise = 1.
    end

    Sshrink = margin_extent(S, ϵ_noise)
    Sexpand = margin_extent(S, -ϵ_noise)
    s_shape = extent_to_shape(Sshrink)
    S_shape = extent_to_shape(Sexpand)

    Q_shape = extent_to_shape(Q_post_extent)
    axes = I+zeros(n_dims, n_dims)

    if ϵ > 0
        # Manual Epsilon Choice
        sd_min = Inf
        for dim_key in keys(Sshrink)
            sd_min = minimum([sd_min, 0.5*abs(Sshrink[dim_key][1]-Sshrink[dim_key][2])])
        end
        shrink_ϵ = minimum([ϵ, sd_min])
        expand_ϵ = ϵ

        Sshrink2 = margin_extent(S, shrink_ϵ)
        Sexpand2 = margin_extent(S, -expand_ϵ)
        s_shape2 = extent_to_shape(Sshrink2)
        S_shape2 = extent_to_shape(Sexpand2)

        shrink_int, _ = does_shape2_contain_shape1(Q_shape, s_shape2, axes=axes)
        expand_int, _ = do_shapes_intersect(Q_shape, S_shape2, axes=axes)
    else
        # Auto Epsilon Choice
        shrink_int, shrink_ϵ = does_shape2_contain_shape1(Q_shape, s_shape, axes=axes)
        expand_int, expand_ϵ = do_shapes_intersect(Q_shape, S_shape, axes=axes)
    end

    if (shrink_int && !expand_int) 
        @info "Shrink Flag: ", shrink_int, "Expand Flag: ", expand_int
        @info "Q Post Extent: ", Q_post_extent 
        @info "Q Post Shape: ", Q_shape
        @info "S Extent: ", S
        @info "S Shrink Ex: ", Sshrink
        @info "S Shrink Shape: ", s_shape
        @info "S Expand Ex: ", Sexpand
        @info "S Expand Shape: " S_shape
        @info "Shrink Ep: ", shrink_ϵ, "Expand Ep: ", expand_ϵ 
        @error "You cannot both interesect and not intersect at the same time!"
    end

    if shrink_int 
        P_bound_lower = calculate_probability_bound_dict(gp_info_dict, σ, shrink_ϵ)
        minPr_LB = P_bound_lower*Pr_noise 
    else
        minPr_LB = 0.
    end

    if !expand_int 
        Pr_bound_upper = calculate_probability_bound_dict(gp_info_dict, σ, expand_ϵ)
        maxPr_UB = 1. - Pr_bound_upper
    else
        maxPr_UB = 1.
    end

    @assert maxPr_UB >= minPr_LB

    if target_region_idx == region_dict.count
        frame_row = [1-maxPr_UB, 1-minPr_LB, [shrink_ϵ, expand_ϵ]]
    else
        frame_row = [minPr_LB, maxPr_UB, [shrink_ϵ, expand_ϵ]]
    end

    return frame_row
end

###############################################################################################
# Local Methods
###############################################################################################
function margin_extent(extent, ep)
    new_dict = Dict()
    for dim_key in keys(extent)
        new_dict[dim_key] = [extent[dim_key][1] + ep, extent[dim_key][2] - ep]
    end
    return return new_dict 
end

function extent_to_shape(extent)
    extent_arrays = [extent[dim_key] for dim_key in keys(extent)]
    vert_tuples = collect(Base.product(extent_arrays...))
    vertices = [hcat([i[j] for j=1:length(i)]...) for i in vert_tuples]
    return Dict("vertices" => vertices[:])
 end

function calculate_probability_bound(gp_info_dict, σ_bounds, ϵ)
    Pr_bound = 1.
    for (i, dim_key) in enumerate(keys(gp_info_dict)) 
        bound_type = gp_info_dict[dim_key].bound_type
        if bound_type == "rkhs-tight" 
            Pr_bound *= 1. - tight_rkhs_bound_prob(σ_bounds[i], ϵ, gp_info_dict[dim_key])

        elseif bound_type == "rkhs-original"
            Pr_bound *= 1. - original_rkhs_bound_prob(σ_bounds[i], ϵ, gp_info_dict[dim_key])
        # elseif bound_type == "munich"
        #     # TODO: Fix this so it works
        #     Pr_bound *= 1. - munich_bound_prob(point, gp_x1, prior_norms[1], ϵ)
        else
            @error "Invalid bound type:" bound_type
        end

    end
    return Pr_bound
end

function calculate_probability_bound_dict(gp_info_dict, σ_bounds, ϵ)
    Pr_bound = 1.
    for dim_key in keys(gp_info_dict) 
        bound_type = gp_info_dict[dim_key].bound_type
        if bound_type == "rkhs-tight" 
            Pr_bound *= 1. - tight_rkhs_bound_prob(σ_bounds[dim_key], ϵ, gp_info_dict[dim_key])
        elseif bound_type == "rkhs-original"
            Pr_bound *= 1. - original_rkhs_bound_prob(σ_bounds[dim_key], ϵ, gp_info_dict[dim_key])
        else
            @error "Invalid bound type:" bound_type
        end
    end
    return Pr_bound
end

function calculate_ϵ_crit(gp_info_dict, σ_bounds)
    candidates = []
    for (i, dim_key) in enumerate(keys(gp_info_dict))
        ϵ = 0.01
        P = 0.0
        bound_type = gp_info_dict[dim_key].bound_type
        while P != 1.0
            if bound_type == "rkhs-tight" 
                P = 1. - tight_rkhs_bound_prob(σ_bounds[i], ϵ, gp_info_dict[dim_key])
            elseif bound_type == "rkhs-original"
                P = 1. - original_rkhs_bound_prob(σ_bounds[i], ϵ, gp_info_dict[dim_key])
            else
                @error "Invalid bound type:" bound_type
            end
            ϵ *= 1.5
        end
        push!(candidates, ϵ)
    end
    return maximum(candidates)
end
