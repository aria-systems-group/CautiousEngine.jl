#= 
Catch all for all the parameters that makes this easy to use;
=#

struct SystemParameters
    mode_tag::String
    unknown_dynamics_fcn
    known_dynamics_fcn
    measurement_noise_dist
    process_noise_dist
    lipschitz_bound::Float64
    dependency_dims::Dict
    n_dims_in::Int
    n_dims_out::Int
end

mutable struct DataParameters 
    data_num::Int
    data_frac_optimize::Float64
    bound_type::String
    noise_sigma::Float64
    epsilon::Float64
    eta::Float64
    safety_dims
    local_radius::Float64
    num_neighbors::Int
end

mutable struct ExperimentParameters 
    experiment_directory::String
    system_tag::String
    specification_file::String
    domain::Dict
    discretization_step::Dict
    verification_mode::String
    random_seed::Int
    system_params::SystemParameters
    data_params
end

function save_experiment_params(params::ExperimentParameters)
	# TODO: Finish this
end

function load_experiment_params(filename::String)
	# TODO: Finish this
end