# Functions to aid in discretizing sets.

" Discretize a given extent into smaller extents with a grid size of delta.
# Arguments
- `set::Dict` - Set to discretize
- `epsilon::Float64` - Size of the discretization
"
function discretize_set(set, grid_sizes)

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
- `epsilon::Float64` - Size of the discretization
"
function discretize_set_lazy(set, epsilon)
    # TODO: Generalize this for n-dimensions
    x = set["x1"]
    y = set["x2"]

    x_epsilons = x[1]:epsilon["x1"]:x[2]
    y_epsilons = y[1]:epsilon["x2"]:y[2]

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
