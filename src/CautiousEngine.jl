# Module that wraps running an experiment into a nice package
module CautiousEngine 

using Logging
using LoggingExtras
using BSON
using Dates
using Distributed
using Base.Threads
using DataFrames
using LinearAlgebra
using Base.Threads
using GaussianProcesses
using Random
using Distributions
using StatsBase
using Serialization
using MAT
using Printf
using IntervalSets
using SparseArrays
using Plots
using Spot

# Lazy sets, the hot new way to work with polytopes
using LazySets

include("metadata.jl")
include("parameters.jl")
include("validation.jl")
include("e2e.jl")

# My own packages
include("GPBounding/GPBounding.jl")
using .GPBounding

include("abstraction/discretization.jl")
include("abstraction/images.jl")
include("abstraction/intersection.jl")
include("abstraction/regions.jl")
include("abstraction/refinement.jl")

include("imdp/dfa.jl")
include("imdp/imdp.jl")
include("imdp/multimodal.jl")
include("imdp/pimdp.jl")
include("imdp/transitions.jl")


include("online/online_control.jl")
include("online/run_simulation.jl")

include("regression/data.jl")
include("regression/gps.jl")
include("regression/error_bounds.jl")

include("synthesis/verification.jl")
include("synthesis/synthesis.jl")

include("visualization/visualize.jl")

export SystemParameters, 
       DataParameters,
       ExperimentParameters

export general_label_fcn

export perform_synthesis_from_result_dirs, 
       generate_transition_bounds, 
       generate_linear_truth, 
       verify_experiment, 
       generate_result_dir_name

## TODO: Where tf to put this
"""
Wrapper function for unpacking the experiment parameters.
"""
function generate_training_data(params::ExperimentParameters)
   
    x_train, y_train = generate_training_data(params.system_params.unknown_dynamics_fcn, params.domain, params.data_params.data_num,
                                              known_fn=params.system_params.known_dynamics_fcn, random_seed=params.random_seed, 
                                              n_dims_out=params.system_params.n_dims_out, process_dist=params.system_params.process_noise_dist, 
                                              measurement_dist=params.system_params.measurement_noise_dist)

    return x_train, y_train
end

end
