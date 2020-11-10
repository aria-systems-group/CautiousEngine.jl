"Calculate the probability that a point is in an extent S given Gaussian processes that model the dynamics.
# Arguments
- `gp_dict::Dict` - A dictionary with Gaussian processes in `:x1` and `:x2` keys
- `point::Float64[]` - A point that the system starts at
- `S::Dict` - The extent of the target set. 
"
function calculate_prob_point_in_S(gp_dict, point, S)
    # Create normal distributions for each component
    normal_dist_x1 = create_normal_dist_from_gp(gp_dict[:x1], point)
    normal_dist_x2 = create_normal_dist_from_gp(gp_dict[:x2], point)
    # Get the intervals of the target set
    Sx_interval = (minimum(S["x1"]), maximum(S["x1"]))
    Sy_interval = (minimum(S["x2"]), maximum(S["x2"]))
    # Calculate the cdf over each region, and multiply to get the total probability
    PrX = calculate_cdf_over_region_1d(normal_dist_x1, Sx_interval)        
    PrY = calculate_cdf_over_region_1d(normal_dist_x2, Sy_interval)        
    Pr = PrX*PrY
    return Pr 
end

"Calculate the probability of being a distance epsilon from a GP at a given point.
# Arguments
- `gp_dict::Dict` - A dictionary with Gaussian processes in `:x1` and `:x2` keys
- `point::Float64[]` - A point that the system starts at
- `epsilon::Float64` - The value of the distance from the process.
"
function calculate_prob_m_epsilon_from_y(gp_dict, point, epsilon)
    # Create normal distributions for each component. Use zero-mean distributions for simpler computation
    normal_dist_x1 = create_normal_dist_from_gp(gp_dict[:x1], point, zero_mean=true)
    normal_dist_x2 = create_normal_dist_from_gp(gp_dict[:x2], point, zero_mean=true)
    # Calculate the CDF over each epsilon interval, and multiply to get the total probability
    PrX = cdf(normal_dist_x1, epsilon) - cdf(normal_dist_x1, -epsilon)
    PrY = cdf(normal_dist_x2, epsilon) - cdf(normal_dist_x2, -epsilon)
    Pr = PrX*PrY
    return Pr 
end

" Create a normal distribution from a Gaussian process and index.
# Arguments
- `gp::GaussianProcess` - The Gaussian process to index
- `point::Float64[]` - A point that the system starts at
"
function create_normal_dist_from_gp(gp, point; zero_mean=false)
    # Get the mean and standard deviation at the index location
    μ, σ = predict_y(gp, hcat(point))
 
    # Create the normal distribution with the std of observations included
    σn = exp(gp.logNoise.value)
    if zero_mean
        dist = Normal(0, sqrt(σ[1]^2 + σn^2))
    else
        dist = Normal(μ[1], sqrt(σ[1]^2 + σn^2))
    end
    return dist
end

" Calculate the cumulative distribution function over a 1d interval.
# Arguments
- `dist::NormalDistribution` - The normal distribution used to calculate the CDF
- `interval::Float64[]` - The extent of the interval to calculate CDF over
"
function calculate_cdf_over_region_1d(dist, interval)
    return cdf(dist, interval[2]) - cdf(dist, interval[1])
end

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

    # i = 0
    # for dim_key in keys(set)
    #     x = set[dim_key]
    #     x_d = x[1]:epsilon:x[2]
    #     num_states = (length(x_d)-1)^length(keys(set))
    #     if i == 0
    #         [[discrete_sets[j+i] = Dict(dim_key =>[x_d[i], x_d[i+1]]) for i=1:length(x_d)-1] for j=0:length(x_d)-1:num_states-1]
    #         i=1
    #     else
    #         [[discrete_sets[j+i][dim_key] = [x_d[i], x_d[i+1]] for i=1:length(x_d)-1] for j=0:length(x_d)-1:num_states-1]
    #     end
    # end
    # x = set["x1"]
    # y = set["x2"]

    # x_epsilons = x[1]:epsilon:x[2]
    # y_epsilons = y[1]:epsilon:y[2]

    
    # idx = 1
    # for i = 1:length(x_epsilons)-1
    #     for j = 1:length(y_epsilons)-1
    #         new_extent = Dict("x1" => [x_epsilons[i], x_epsilons[i+1]], "x2" => [y_epsilons[j], y_epsilons[j+1]])
    #         discrete_sets[idx] = new_extent
    #         idx += 1
    #     end 
    # end
    return discrete_sets
end

" Discretize a given extent into smaller extents with a grid size of delta.
# Arguments
- `set::Dict` - Set to discretize
- `epsilon::Float64` - Size of the discretization
"
function discretize_set_new(set, epsilon)
    x = set["x1"]
    y = set["x2"]

    x_epsilons = x[1]:epsilon:x[2]
    y_epsilons = y[1]:epsilon:y[2]

    discrete_sets = Dict()
    discrete_extents = Dict()
    idx = 1
    for i = 1:length(x_epsilons)-1
        for j = 1:length(y_epsilons)-1
            new_extent = Extent("x1" => [x_epsilons[i], x_epsilons[i+1]], "x2" => [y_epsilons[j], y_epsilons[j+1]])
            new_polygon = VPolygon([[x_epsilons[i], y_epsilons[j]], [x_epsilons[i], y_epsilons[j+1]], 
                                    [x_epsilons[i+1], y_epsilons[j+1]], [x_epsilons[i+1], y_epsilons[j]]])
            discrete_sets[idx] = new_polygon 
            discrete_extents[idx] = new_extent
            idx += 1
        end 
    end
    return discrete_sets, discrete_extents
end

" Expands or shrinks a rectangle according to parameter epsilon.
# Arguments
- `rect::Dict` - Dictionary containing shape to expand or shrink
- `epsilon::Float64` - Size of the expansion or shrinkage
"
function margin_rectangle(rect, epsilon)
    x = rect["x1"]
    y = rect["x2"]

    new_x_min = minimum(x) + epsilon   
    new_x_max = maximum(x) - epsilon
    new_y_min = minimum(y) + epsilon
    new_y_max = maximum(y) - epsilon

    # TODO: Remove the redunant points
    margined_rectangle = Dict("x1" => [new_x_min, new_x_min, new_x_max, new_x_max],
                              "x2" => [new_y_min, new_y_max, new_y_max, new_y_min])
    return margined_rectangle
end
