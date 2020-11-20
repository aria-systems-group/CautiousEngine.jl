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
