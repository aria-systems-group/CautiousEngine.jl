using Random
using Distributions

function generate_data_bounded_gaussian(params, dyn_fn)
    domain = params.domain
    m = params.data_params.data_num
    noise_sigma = params.data_params.noise_sigma
    n_dims_in = length(keys(domain))
    n_dims_out = length(dyn_fn(zeros(n_dims_in, 1)))
    # x1_train = rand(Uniform(domain["x1"][1], domain["x1"][2]), m)
    # x2_train = rand(Uniform(domain["x2"][1], domain["x2"][2]), m)
    mt = MersenneTwister(params.random_seed)
    x_train = hcat([rand(mt, Uniform(domain["x$i"][1], domain["x$i"][2]), m) for i=1:n_dims_in]...)
    y_train = mapslices(dyn_fn, x_train, dims=2) + rand(mt, Truncated(Normal(0, noise_sigma), -noise_sigma, noise_sigma), (m,n_dims_out)) 
    return x_train, y_train
end

function generate_data_gaussian(X, dyn_fn, m, noise_sigma)
    x1_train = rand(Uniform(X["x1"][1], X["x1"][2]), m)
    x2_train = rand(Uniform(X["x2"][1], X["x2"][2]), m)
    x_train = [x1_train x2_train]
    y_train = mapslices(dyn_fn, x_train, dims=2) + rand(Normal(0, noise_sigma), (m,2)) 
    return x_train, y_train
end
