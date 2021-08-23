# Contains functions for generating synthetic training data

"""
Generates synthetic data from an unknown function by uniformly sampling datapoints over the domain. Returns the input and output points. 

    Note that you cannot use the `measurement_dist` and `known_fn` together.

# Arguments 
- `unknown_fn`: the handle of the unknown function.
- `domain::Dict`: the dictionary indicating the state space domain with keys `"x1","x2"...` and corresponding min-max bounds.
- `data_num::Int`: the number of datapoints to generate.
- `known_fn=nothing`: the known part of the function.
- `random_seed::Int=11`: random seed for the uniform sample (and noise processes, if any).
- `n_dims_out::Int=-1`: dimension of the output of the function, defaults to same dimension as input.
- `process_dist=nothing`: process noise distristribution.
- `measurement_dist=nothing`: measurement noise distristribution.
"""
function generate_training_data(unknown_fn, domain::Dict, data_num::Int; 
                                known_fn=nothing, random_seed::Int=11, n_dims_out::Int=-1, process_dist=nothing, measurement_dist=nothing)

    if !isnothing(known_fn)
        f_sub = (x) -> unknown_fn(x) + known_fn(x)
    else 
        f_sub = (x) -> unknown_fn(x)
    end

    n_dims_in = length(domain)
    n_dims_out = n_dims_out > 0 ? n_dims_out : n_dims_in 
    mt = MersenneTwister(random_seed)

    x_train = vcat([rand(mt, Uniform(domain["x$i"][1], domain["x$i"][2]), 1, data_num) for i=1:n_dims_in]...)
    y_train = mapslices(f_sub, x_train, dims=1) 
   
    # Account for process noise...
    if !isnothing(process_dist)
        y_train += rand(mt, process_dist, (n_dims_out, data_num)) 
    end

    # Account for measurement noise...
    if !isnothing(measurement_dist)
        y_train += rand(mt, measurement_dist, (n_dims_out, data_num)) 
    end

    # Account for the known part, if any...
    if !isnothing(known_fn)
        y_train = y_train - mapslices(known_fn, x_train, dims=1)
    end

    @assert size(x_train, 2) == size(y_train, 2)

    return x_train, y_train
end

function save_training_data(x_train, y_train, exp_dir::String)
    training_data = Dict()
    training_data[:input] = x_train
    training_data[:output] = y_train
    bson("$exp_dir/training_data.bson", training_data)
end

function load_training_data(filename::String)
    training_data = BSON.load(filename)
    return training_data[:input], training_data[:output]
end

