# Functions to aid in discretizing sets.

" Discretize a given extent into smaller extents with a grid size of delta.
# Arguments
- `set::Dict` - Set to discretize
- `grid_sizes::Dict` - Discretization delta for each component
"
function discretize_set(set::Dict, grid_sizes::Dict)

    extents = []
    for dim_key in keys(set)
        x = set[dim_key]
        x_d = x[1]:grid_sizes[dim_key]:x[2]
        push!(extents, [[x_d[i], x_d[i+1]] for i=1:length(x_d)-1])
    end

    state_extents = Base.product(extents...)
    discrete_sets = Dict()
    for (i, state) in enumerate(state_extents)
        for (j, dim_key) in enumerate(collect(keys(set)))
            if j==1
                discrete_sets[i] = Dict(dim_key => state[j]) 
            else
                discrete_sets[i][dim_key] = state[j]
            end
        end
    end

    return discrete_sets
end

" Discretize a given extent into smaller extents with a grid size of delta.
# Arguments
- `set::Dict` - Set to discretize
- `grid_sizes::Dict` - Discretization delta for each component
"
function discretize_set_lazy(set::Dict, grid_sizes::Dict)
    # TODO: Generalize this for n-dimensions
    x = set["x1"]
    y = set["x2"]

    x_epsilons = x[1]:grid_sizes["x1"]:x[2]
    y_epsilons = y[1]:grid_sizes["x2"]:y[2]

    discrete_sets = Dict()
    discrete_extents = Dict()
    idx = 1
    for i = 1:length(x_epsilons)-1
        for j = 1:length(y_epsilons)-1
            new_extent = Dict("x1" => [x_epsilons[i], x_epsilons[i+1]], "x2" => [y_epsilons[j], y_epsilons[j+1]])
            new_polygon = VPolygon([[x_epsilons[i], y_epsilons[j]], [x_epsilons[i], y_epsilons[j+1]], 
                                    [x_epsilons[i+1], y_epsilons[j+1]], [x_epsilons[i+1], y_epsilons[j]]])
            discrete_sets[idx] = new_polygon 
            discrete_extents[idx] = new_extent
            idx += 1
        end 
    end
    return discrete_sets, discrete_extents
end

" Refine a discrete state uniformly and update the region info dict.
# Arguments
- `region_id::Int` - Region id to refine 
- `region_info_dict:Dict` - Dictionary containing all of the region info
"
function uniform_refinement(region_id, extent_dict)
    extent = extent_dict[region_id]
    delta_dict = Dict()
    for dim in keys(extent)
        delta_dict[dim] = (extent[dim][2] - extent[dim][1])/2
    end

    refined_sets = discretize_set(extent_dict[region_id], delta_dict)

    # ! Does ordering matter? I don't think it should
    N_old = maximum(keys(extent_dict))
    last_extent = extent_dict[N_old]

    # Remove old state
    pop!(extent_dict, region_id)
    new_keys = []
    for i=1:refined_sets.count
        new_key = N_old + i - 1 
        extent_dict[N_old + i - 1] = refined_sets[i]
        new_keys = new_keys ∪ [new_key]
    end

    extent_dict[N_old + refined_sets.count + 1] = last_extent

    return extent_dict, new_keys
end