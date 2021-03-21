# Currently a placeholder that will hold data-generation related code.

function generate_training_data(params)
    if !isnothing(params.system_params.known_dynamics_fcn)
        f_sub = (x) -> params.system_params.unknown_dynamics_fcn(x) + params.system_params.known_dynamics_fcn(x)
    else 
        f_sub = (x) -> params.system_params.unknown_dynamics_fcn(x)
    end

    domain = params.domain
    m = params.data_params.data_num
    noise_sigma = params.data_params.noise_sigma
    n_dims_in = params.system_params.n_dims_in 
    n_dims_out = params.system_params.n_dims_out 
    mt = MersenneTwister(params.random_seed)

    # x_train = hcat([rand(mt, Uniform(domain["x$i"][1], domain["x$i"][2]), m) for i=1:n_dims_in]...)
    x_train = vcat([rand(mt, Uniform(domain["x$i"][1], domain["x$i"][2]), 1, m) for i=1:n_dims_in]...)
    @info x_train
    y_train = mapslices(f_sub, x_train, dims=1) 
   
    # Account for process noise...
    if !isnothing(params.system_params.process_noise_dist)
        # y_train += rand(mt, params.system_params.process_noise_dist, (m,n_dims_out)) 
        y_train += rand(mt, params.system_params.process_noise_dist, (n_dims_out, m)) 
    end

    # Account for measurement noise...
    if !isnothing(params.system_params.measurement_noise_dist)
        y_train += rand(mt, params.system_params.measurement_noise_dist, (n_dims_out, m)) 
    end

    # Account for the known part, if any...
    # TODO: I don't think we actually need this part.
    if !isnothing(params.system_params.known_dynamics_fcn)
        y_train = y_train - mapslices(params.system_params.known_dynamics_fcn, x_train, dims=1)
    end

    @assert size(x_train, 2) == size(y_train, 2)

    return x_train, y_train
end