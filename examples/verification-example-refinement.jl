# Include the known function
using Random
using Distributions
using IntervalSets
using SparseArrays

using CautiousEngine

# Get the absolute path of the experiment dir
EXPERIMENT_DIR = @__DIR__ 
experiment_type = "verification-example-refinement" 
global_exp_dir = "$EXPERIMENT_DIR/$experiment_type/global"
specification_file = "foo"

#======================================================================================================
0. Define the Unknown Function, Process Noise, and the Full Dynamics Function
======================================================================================================#
system_mat = [0.9 0.5; -0.5 0.4]
f = (x) -> system_mat*x 

system = "mode1"
known_part = nothing        

# Define the zero-mean process noise
σ_proc = 0.01
measurement_noise_dist = Normal(0., σ_proc)             
process_noise_dist = nothing                                # Process noise is only used in synthesis, not verification 

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
number_of_datapoints = 400                          # Number of data points per mode to generate
run_exps_flag = true                                # Set to true to rerun the single-mode constructions 
m_opt = -1                                          # Unused
bound_type = "rkhs-tight"                           # GP bound type 
epsilon = -1.                                       # Distance bound parameter, set to -1 for automatic determination

Lf_bound = 1.5                                 # Lipschitz constant bound for the mode over X
random_seed = 11                             	# Random seed for the random twister 
eta = -1.                                       # η parameter set in a naiive manner

# Setup experiment parameter structures (defined in CautiousEngine.jl)
system_params = SystemParameters(system, f, known_part, measurement_noise_dist, process_noise_dist, Lf_bound, dependency_dims, 2, 2)
data_params = DataParameters(number_of_datapoints, m_opt, bound_type, σ_proc, epsilon, eta, safety_dims, -1., -1.)
experiment_params = ExperimentParameters(exp_dir, experiment_type, specification_file, X, grid_sizes, "foo", random_seed, system_params, data_params)

_, trans_mats = CautiousEngine.end_to_end_transition_bounds(experiment_params, reuse_gps_flag=true, reuse_regions_flag=true)
results_path = experiment_params.experiment_directory

#======================================================================================================
2.  Construct the IMDP Labels and call the appropriate verification solver wrapper. 
======================================================================================================#
# Wrapper function
imdp = CautiousEngine.create_simple_imdp(trans_mats["minPr"], trans_mats["maxPr"]) 
refinement_steps = 3

k_vals = [1,]
for k in k_vals
    global imdp
    verification_result_mat = CautiousEngine.Globally(imdp, "safe", k, "$results_path/imdp.txt")
    CautiousEngine.save_legacy_mats(verification_result_mat, results_path, k)
    CautiousEngine.plot_2d_verification_results(results_path, prob_plots=true)
    working_dir = results_path

    # TODO: Reusing regions only works for a single instance of k. Handle this in a smarter way.
    refine_result_dirs = CautiousEngine.safety_based_refinement(experiment_params, results_path, verification_result_mat, refinement_steps, horizon=k, reuse_regions_flag=true, reuse_transition_mats_flag=false)

    for res in refine_result_dirs
        CautiousEngine.plot_2d_verification_results(res, prob_plots=true, state_outlines_flag=false)
    end
end
