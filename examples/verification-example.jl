using Distributions
using Random
using SparseArrays 

using CautiousEngine

# Get the absolute path of the experiment dir
EXPERIMENT_DIR = @__DIR__ 
experiment_type = "verification-example" 
global_exp_dir = "$EXPERIMENT_DIR/$experiment_type/global"
specification_file = "foo"

#======================================================================================================
0. Define the Unknown Function, Process Noise, and the Full Dynamics Function
======================================================================================================#
# Full dynamics function with three control modes
A1 = [0.8 0.5; 0. 0.5]
A2 = [0.5 0.; -0.5 0.8]
g1 = (x) -> A1*x
g2 = (x) -> A2*x
g3 = (x) -> [x[1]+(sin(x[2])); 0.85*sqrt(abs(x[1]))*x[2]]
unknown_mode_list = [g1, g2, g3]
mode_names = ["mode1", "mode2", "mode3"]
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
grid_sizes = Dict("x1" => 0.25, "x2" => 0.25)       # Specifies the discretization coarseness in each dimension
safety_dims = ["x1", "x2"]                          # Dimensions that have safety bounds
dependency_dims = Dict("x1" => [1, 1],              # Defines dependencies on other states  
                       "x2" => [1, 1])
number_of_datapoints = 100                          # Number of data points per mode to generate
run_exps_flag = true                                # Set to true to rerun the single-mode constructions 
m_opt = -1                                          # Unused
bound_type = "rkhs-tight"                           # GP bound type 
epsilon = -1.                                       # Distance bound parameter, set to -1 for automatic determination

for (i, mode) in enumerate(unknown_mode_list)
    system = mode_names[i]                          # Get the mode name
    Lf_bound = 1.0                                  # Lipschitz constant bound for the mode over X
    random_seed = 11 + i                            # Random seed for the random twister 
    eta = σ_proc                                    # η parameter set in a naiive manner

    @info system, number_of_datapoints, σ_proc, bound_type, Lf_bound, random_seed 

    # Setup experiment parameter structures (defined in CautiousEngine.jl)
    system_params = SystemParameters(system, mode, known_part, measurement_noise_dist, process_noise_dist, Lf_bound, dependency_dims, 2, 2)
    data_params = DataParameters(number_of_datapoints, m_opt, bound_type, σ_proc, epsilon, eta, safety_dims, -1., -1.)
    experiment_params = ExperimentParameters(exp_dir, experiment_type, specification_file, X, grid_sizes, "foo", random_seed, system_params, data_params)

    # If constructing the single-mode systems, run it! 
    if run_exps_flag
        CautiousEngine.end_to_end_transition_bounds(experiment_params, single_mode_verification=true, reuse_regions_flag=false)
    end

    # Construct the directory where results will be
    res_dir = CautiousEngine.create_experiment_directory(experiment_params)
    @info "Mode results can be found in $res_dir" 
    push!(res_dirs, res_dir)
end

#======================================================================================================
2.  Construct the IMDP Labels and call the appropriate verification solver wrapper. 
======================================================================================================#
default_label = "safe"                      # Default label of states
unsafe_label = "!safe"                      # The label of unsafe states
labels_dict = Dict("safe" => [X])           # Dictionary with key if label and value of corresponding extents 
# Wrapper function
lbl = (state_extent; unsafe=false) -> general_label_fcn(state_extent, default_label, unsafe_label, labels_dict, unsafe=unsafe)
switched_system_dir = "$exp_dir/m$number_of_datapoints" 
imdp = CautiousEngine.create_imdp_from_result_dirs(res_dirs, switched_system_dir, label_fn=lbl)

k_vals = [1, 2, 3, 4, 5, -1]                # Verification horizon (-1 means unbounded)
for k in k_vals
    verification_result_mat = CautiousEngine.Globally(imdp, "safe", k, "$switched_system_dir/switched-system/imdp-$k.txt")
    CautiousEngine.save_legacy_mats(verification_result_mat, "$switched_system_dir/switched-system", k)
end
CautiousEngine.plot_2d_verification_results("$switched_system_dir/switched-system")
