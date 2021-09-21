using Distributions
using Random
using IntervalSets
using SparseArrays

using CautiousEngine

# Get the absolute path of the experiment dir
EXPERIMENT_DIR = @__DIR__ 
experiment_type = "safe+live_example" 
global_exp_dir = "$EXPERIMENT_DIR/$experiment_type/global"
online_exp_dir = "$EXPERIMENT_DIR/$experiment_type/online"
specification_file = "reach-specification.jl"

#======================================================================================================
0. Define the Unknown Function, Process Noise, and the Full Dynamics Function
======================================================================================================#
# Four control modes defined - g is the unknown function, f is the total (including known) dynamics fcn
v = 0.5
g1 = (x) -> [v + 0.1*sin(x[2]), 0.2*cos(x[1])]*0.5
g2 = (x) -> [-v + 0.1*sin(x[2]), 0.2*cos(x[1])]*0.5
g3 = (x) -> [0.2*cos(x[2]), v + 0.1*sin(x[1])]*0.5
g4 = (x) -> [0.2*cos(x[2]), -v + 0.1*sin(x[1])]*0.5

# Define the zero-mean process noise
σ_proc = 0.01
process_noise_dist = Truncated(Normal(0., σ_proc), -σ_proc, σ_proc)
measurement_noise_dist = nothing                    # Measurement noise is only used in verification, not synthesis

# Full dynamics function with four control modes
f1 = (x) -> x + g1(x) + rand(process_noise_dist, (2,1))
f2 = (x) -> x + g2(x) + rand(process_noise_dist, (2,1)) 
f3 = (x) -> x + g3(x) + rand(process_noise_dist, (2,1))
f4 = (x) -> x + g4(x) + rand(process_noise_dist, (2,1))

mode_list = [f1, f2, f3, f4]
mode_names = ["mode1", "mode2", "mode3", "mode4"]
unknown_mode_list = [g1, g2, g3, g4]
known_part = (x) -> [x[1], x[2]]                    # Known part of the dynamics function, can only be this currently.

#======================================================================================================
1.  Perform GP regression and construct IMDP for each unknown mode.
======================================================================================================#
exp_dir = global_exp_dir
res_dirs = Array{String,1}()                        # Array that holds the single-mode result directories
X = Dict("x1" => [-2., 2.], "x2" => [-2., 2.])      # Defines the compact safe set
grid_sizes = Dict("x1" => 1.0, "x2" => 1.0)       # Specifies the discretization coarseness in each dimension
safety_dims = ["x1", "x2"]                          # Dimensions that have safety bounds
dependency_dims = Dict("x1" => [1, 1],              # Defines dependencies on other states  
                       "x2" => [1, 1])
number_of_datapoints = 500                           # Number of data points per mode to generate
run_exps_flag = true                                # Set to true to rerun the single-mode constructions 
m_opt = -1                                          # Unused
bound_type = "rkhs-tight"                           # GP bound type 
epsilon = -1.                                       # Distance bound parameter, set to -1 for automatic determination

params = []
for (i, mode) in enumerate(unknown_mode_list)
    system = mode_names[i]                          # Get the mode name
    Lf_bound = 1.0                                  # Lipschitz constant bound for the mode over X
    random_seed = 11 + i                            # Random seed for the random twister 
    eta = σ_proc                                    # η parameter set in a naiive manner

    @info system, number_of_datapoints, σ_proc, bound_type, Lf_bound, random_seed 

    # Setup experiment parameter structures (defined in CautiousEngine.jl)
    system_params = SystemParameters(system, mode, known_part, measurement_noise_dist, process_noise_dist, Lf_bound, dependency_dims, 2, 2)
    data_params = DataParameters(number_of_datapoints, m_opt, bound_type, σ_proc, epsilon, eta, safety_dims, -1., -1.)
    experiment_params = ExperimentParameters(exp_dir, experiment_type, specification_file, X, grid_sizes, specification_file, random_seed, system_params, data_params)
    push!(params, experiment_params)
    # If constructing the single-mode systems, run it! 
    if run_exps_flag
        CautiousEngine.end_to_end_transition_bounds(experiment_params, reuse_gps_flag=false, reuse_regions_flag=false)
    end

    # Construct the directory where results will be
    res_dir = experiment_params.experiment_directory
    @info "Mode results can be found in $res_dir" 
    push!(res_dirs, res_dir)
end

#======================================================================================================
2.  Construct the safety strategy
======================================================================================================#
# if false
    refinement_steps = 4 
    synth_dir = exp_dir                                 # Specify where the synthesis result should be saved 
    !isdir(synth_dir) && mkdir(synth_dir)        

    safe_extents = [Dict("x1" => -2.0.. 2.0,             # List containing extents corresponding to a label 
                        "x2" => -2.0.. 2.0)]
    labels_dict = Dict("safe" => safe_extents)           # Dictionary with key if label and value of corresponding extents 
    default_label = "safe"
    unsafe_label = "!safe"

    # Wrapper function
    lbl = (state_extent; unsafe=false) -> general_label_fcn(state_extent, default_label, unsafe_label, labels_dict, unsafe=unsafe)
    system_tag = "synthesis-example-m$number_of_datapoints" # Tag for this example
    # TODO: Change horizon
    horizon = 1 
    refinement_result_dirs = CautiousEngine.safety_based_synthesis(params, res_dirs, "$synth_dir/$system_tag", system_tag, lbl, refinement_steps, horizon=horizon)

    for res in refinement_result_dirs
        CautiousEngine.plot_2d_verification_results(res, prob_plots=true, state_outlines_flag=true)
    end
# end

#======================================================================================================
3.  Construct the liveness strategy
======================================================================================================#
refinement_steps = 3
synth_dir = exp_dir                                 # Specify where the synthesis result should be saved 
!isdir(synth_dir) && mkdir(synth_dir)        

des_extents = [Dict("x1" => -1.0.. 1.0,             # List containing extents corresponding to a label 
                    "x2" => -1.0.. 1.0)]
labels_dict = Dict("a" => des_extents)           # Dictionary with key if label and value of corresponding extents 
default_label = "!a"
# TODO: Don't need unsafe label!
unsafe_label = "!a"

# Wrapper function
lbl = (state_extent; unsafe=false) -> general_label_fcn(state_extent, default_label, unsafe_label, labels_dict, unsafe=unsafe)
# TODO: Change horizon
horizon = -1
refinement_result_dirs = CautiousEngine.spec_based_refinement(params, res_dirs, "$synth_dir", system_tag, lbl, refinement_steps, horizon=horizon) 

for res in refinement_result_dirs
    CautiousEngine.plot_2d_verification_results(res, prob_plots=true, state_outlines_flag=true, num_dfa_states=2)
end