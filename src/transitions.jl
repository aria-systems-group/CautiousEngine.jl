"Calculate the transition probabilities between a pair of regions."
function region_pair_transitions(region_pair, region_dict, region_post_dict, ϵ, gp_info_dict, params)
    dict_row = region_post_dict[region_pair[1]]
    σ_bounds = dict_row[2] 

    # Get the subset of the safe region that matters.
    safe_keys = params.data_params.safety_dims
    # if length(safe_keys) == 0 && region_pair[2] == -11
    #     # Everything is safe
    #     return [region_pair[1], region_pair[2],  1.0, 1.0, [-1., -1.], [0., 0.]] 
    # end
    S = Dict()
    Q_post_extent = Dict()
    for dim_key in safe_keys
        S[dim_key] = region_dict[region_pair[2]][dim_key] 
        Q_post_extent[dim_key] = dict_row[1][dim_key]
    end

    # For now, do a naiive shrink - first noise, then the next
    # TODO: Pass in real noise value
    noise_dist = TruncatedNormal(0, params.data_params.noise_sigma, -0.01, 0.01)
    ϵ_noise = params.data_params.eta
    Pr_noise = (cdf(noise_dist, ϵ_noise) - cdf(noise_dist, -ϵ_noise))^length(σ_bounds)

    function margin_extent(extent, ep)
        new_dict = Dict()
        for dim_key in keys(extent)
            new_dict[dim_key] = [extent[dim_key][1] + ep, extent[dim_key][2] - ep]
        end
        return return new_dict 
    end

    Sshrink = margin_extent(S, ϵ_noise)
    Sexpand = margin_extent(S, -ϵ_noise)

    s_shape = extent_to_shape(Sshrink)
    S_shape = extent_to_shape(Sexpand)
    Q_shape = extent_to_shape(Q_post_extent)
    n_dims = length(safe_keys)
    axes = I+zeros(n_dims, n_dims)

    if ϵ < 0
        # ϵ_lower = lower_bound_ϵ(Q_post_extent, Sshrink, abs(S["x1"][1] - S["x1"][2])) 
        # ϵ_upper = upper_bound_ϵ(Q_post_extent, Sexpand)
        shrink_int, shrink_ϵ = does_shape2_contain_shape1(Q_shape, s_shape, axes=axes)
        expand_int, expand_ϵ = do_shapes_intersect(Q_shape, S_shape, axes=axes)
        # @assert !(shrink_int && !expand_int)
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
    else
        shrink_ϵ = ϵ
        expand_ϵ = ϵ
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
function extent_to_shape(extent)
    extent_arrays = [extent[dim_key] for dim_key in keys(extent)]
    vert_tuples = collect(Base.product(extent_arrays...))
    vertices = [hcat([i[j] for j=1:length(i)]...) for i in vert_tuples]
    return Dict("vertices" => vertices[:])
 end

# function lower_bound_ϵ(Q, S, ϵ_max)
#     # If all of these checks are false, then Q is necessarily inside of S and !(...) is true.
#     I_S = !(Q["x1"][1] < S["x1"][1] || Q["x1"][2] > S["x1"][2] || Q["x2"][1] < S["x2"][1] || Q["x2"][2] > S["x2"][2]) #and_point_of_Q_post_outside_of_S_small
#     if I_S
#         # How much can we shrink S until there is a point in Q that is outside...
#         ϵ = minimum([abs(Q["x1"][1] - S["x1"][1]), abs(Q["x1"][2] - S["x1"][2]), abs(Q["x2"][1] - S["x2"][1]), abs(Q["x2"][2] - S["x2"][2]), ϵ_max]) 
#     else
#         # If Q has a point outside of S, epsilon is arbitrary. 
#         ϵ = -1.
#     end
#     return ϵ
# end

# function upper_bound_ϵ(Q, S)
#     # If all of these checks are true, then Q necessarily intersects S and (...) is true.
#     I_S = (Q["x1"][1] < S["x1"][2] && Q["x1"][2] > S["x1"][1] && Q["x2"][2] > S["x2"][1] && Q["x2"][1] < S["x2"][2])
#     if I_S
#         # Epsilon is arbitrary
#         ϵ = -1.
#     else
#         # How much can we expand S until we get an intersection?
#         y_min = minimum([abs(Q["x2"][1] - S["x2"][1]), abs(Q["x2"][2] - S["x2"][2]), abs(Q["x2"][2] - S["x2"][1]), abs(Q["x2"][1] - S["x2"][2])]) 
#         x_min = minimum([abs(Q["x1"][1] - S["x1"][1]), abs(Q["x1"][2] - S["x1"][2]), abs(Q["x1"][2] - S["x1"][1]), abs(Q["x1"][1] - S["x1"][2])]) 
#         ϵ = maximum([y_min, x_min])
#         @assert ϵ > 0
#     end
#     return ϵ
# end

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
    # @info Pr_bound, ϵ
    return Pr_bound
end
