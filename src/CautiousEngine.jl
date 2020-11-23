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
using MATLAB
using MAT
using Printf
using IntervalSets
using Plots

# Lazy sets, the hot new way to work with polytopes
using LazySets

# My own packages
using GPBounding
include("transitions.jl")
include("error_bounds.jl")
include("discretization.jl")
include("data.jl")
include("verification.jl")
include("create_dot_graphs.jl")
include("imdp.jl")
include("dfa.jl")
include("pimdp.jl")
include("multimodal.jl")
include("visualize.jl")
include("intersection.jl")
include("run_simulation.jl")

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
    data_num_optimize::Int
    bound_type::String
    noise_sigma::Float64
    epsilon::Float64
    eta::Float64
    safety_dims
end

# TODO: Switch to params structure
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

struct GPInfo
    gp
    γ_bound::Float64
    RKHS_bound::Float64
    bound_type::String
    post_scale_factor::Float64
    Kinv
end

export SystemParameters, DataParameters, ExperimentParameters

export perform_synthesis_from_result_dirs, generate_transition_bounds, generate_linear_truth, verify_experiment, generate_result_dir_name

export plot_results_from_file, plot_gp_field_slice

# TODO: Fix this function
function generate_linear_truth(params, system_matrix; single_mode_verification=false)
    logfile = initialize_log(params)
    timing_info = Dict()
    total_runtime = 0.
    # TODO: Move this
    @info "Determining the post-images of the regions under the linear map."
    region_time = @elapsed begin 
        region_dict, region_pairs, extents_dict = create_region_data_new(params.domain, params.discretization_step)

        # Here, we can calculate the mean and covariance at lots of sampled points in each region
        region_post_dict = Dict()

        L = length(keys(region_dict))
        for i=1:length(keys(region_dict)) 
            # TODO: Fix the indeces to not use -11
            if i == L 
                region = region_dict[-11]
            else
                region = region_dict[i]
            end
            
            # Lazy representation saves time
            # TODO: Add the process noise term here for verification
            region_post = LinearMap(0.9999*system_matrix, region) 

            if i == L 
                region_post_dict[-11] = region_post 
            else
                region_post_dict[i] = region_post 
            end
        end 

        region_data = Dict()
        region_data[:pairs] = region_pairs
        region_data[:extents] = extents_dict 
        region_data[:posts] = region_post_dict 
    end
    timing_info["region_bound_time_s"] = region_time 
    total_runtime += region_time
    @info "Region generation time: " region_time
    save_region_data(params, region_data)
    @info "Generating the transition bounds..."
    bound_time = @elapsed begin
        @info "Calculating transition probability bounds between regions..."
        results_df = DataFrame(Set1 = Int[], Set2 = Int[], MinPr = Float64[], MaxPr = Float64[], MinPrPoint = Array[], MaxPrPoint = Array[])
        for region_pair in region_pairs
            r1 = region_dict[region_pair[1]]
            r2 = region_dict[region_pair[2]]
            r1_post = region_post_dict[region_pair[1]]
        
            # Check lazy intersection
            cap = Intersection(r1_post, r2)
            if isempty(cap)
                prange = [0., 0.]
            else
                if isequivalent(r1_post, cap)
                    prange = [1., 1.]
                else
                    prange = [0., 1.]
                end
            end   

            df_row = [region_pair[1], region_pair[2], prange[1], prange[2], [-1.], [-1.]]
            push!(results_df, df_row)
        end
        num_states = length(keys(region_dict))
        minPr_mat = zeros((num_states, num_states))
        maxPr_mat = zeros((num_states, num_states))
        for i = 1:1:num_states
            sub_row = results_df[results_df.Set1 .== i, :]
            if i == num_states
                minPr_mat[i, i] = 1.
                maxPr_mat[i, i] = 1.
            else
                for j = 1:1:num_states
                    if j == num_states
                        subsub_row = sub_row[sub_row.Set2 .== -11, :]
                        minPr_mat[i, j] = 1. - subsub_row.MaxPr[1]
                        maxPr_mat[i, j] = 1. - subsub_row.MinPr[1]
                    else
                        subsub_row = sub_row[sub_row.Set2 .== j, :]
                        minPr_mat[i, j] = subsub_row.MinPr[1]
                        maxPr_mat[i, j] = subsub_row.MaxPr[1]
                    end
                end
            end
        end

        result_mats = Dict("minPr" => minPr_mat, "maxPr" => maxPr_mat)
    end
    total_runtime += bound_time
    timing_info["transition_bound_time_s"] = bound_time 
    @info "Bound generation time: " bound_time
    save_transition_matrices(params, result_mats)

    if single_mode_verification
        @info "Performing safety verification on single mode..."
        verification_time = @elapsed begin 
            # script_path = "/Users/john/Projects/Julia/BMDPVerification"
            # run_MATLAB_verification(script_path, create_experiment_directory(params), -1, params.verification_mode)
            perform_imdp_verification(params, result_mats, 1)
        end
        timing_info["single_verification_time_s"] = verification_time  
        @info "Verification time: " verification_time
    end

    @info "Total runtime: " total_runtime
    timing_info["total_runtime_s"] = total_runtime
    save_time_info(params, timing_info)
    
    flush(logfile)
    close(logfile)

    return timing_info

    # Holder of the metadata rows
    

    # TODO: Save time info to log file
    

    # Save the matrices for use in MATLAB verification tool
    save_mats = false
    if save_mats
        @info "Saving matrices in MATLAB format..."
        num_states = length(keys(region_dict))
        A_min = zeros((num_states, num_states))
        A_max = zeros((num_states, num_states))
        for i = 1:1:num_states
            sub_row = results_df[results_df.Set1 .== i, :]
            if i == num_states
                A_min[i, num_states] = 1.
                A_max[i, num_states] = 1.
            else
                for j = 1:1:num_states
                    if j == num_states
                        subsub_row = sub_row[sub_row.Set2 .== -11, :]
                        A_min[i, j] = 1. - subsub_row.MaxPr[1]
                        A_max[i, j] = 1. - subsub_row.MinPr[1]
                    else
                        subsub_row = sub_row[sub_row.Set2 .== j, :]
                        A_min[i, j] = subsub_row.MinPr[1]
                        A_max[i, j] = subsub_row.MaxPr[1]
                    end
                end
            end
        end

        matwrite("$exp_dir/transition_mats.mat", Dict("minPr" => A_min, "maxPr" => A_max))

    end
  
end

# TODO: This will eventually be the end-to-end function
# function perform_synthesis_experiment(params, dynamics_fcn)
   
#     generate_transition_bounds(params, dynamics_fcn)
#     # Now do all the synthesis stuff here.
#     perform_synthesis_from_result_dirs(res_dirs, params.experiment_directory, labels, params.system_tag, params.specification_file) 

# end

function perform_synthesis_from_result_dirs(res_dirs, exp_dir, system_tag, spec_file, label_fcn; 
                                            plot_graphs=true, rerun_flag=false, modes=nothing, add_opt=false,
                                            region_specific="", num_sims=10, sim_points=nothing)

    @info "Creating the IMDP..."
    imdp = create_imdp_from_result_dirs(res_dirs, "$exp_dir/$system_tag")    
    regions_file = @sprintf("%s/regions%s.bson", res_dirs[1], region_specific)
    create_imdp_labels(label_fcn, imdp, regions_file)
    # imdp.labels[length(keys(imdp.labels))] = "!a∧b"
    plot_graphs ? create_graph_from_imdp(imdp, "$exp_dir/$system_tag/imdp.gv") : nothing

    @info "Creating the DFA..."
    spec_tag = basename(spec_file)[1:end-3]
    include("$spec_file")
    dfa = DFA(dfa_states, dfa_props, dfa_transitions, dfa_accepting_state, dfa_sink_state, dfa_initial_state)
    dst_dir = "$exp_dir/$system_tag/$spec_tag"
    isdir(dst_dir) ? nothing : mkpath(dst_dir)
    plot_graphs ? create_graph_from_dfa(dfa, "$dst_dir/specification.gv") : nothing

    @info "Creating pimdp..."
    pimdp = construct_DFA_IMDP_product(dfa, imdp)
    plot_graphs ? create_graph_from_pimdp(pimdp, "$dst_dir/pimdp.gv") : nothing

    pimdp_filename = "$dst_dir/pimdp.txt"
    write_pimdp_to_file(pimdp, pimdp_filename)

    @info "Synthesizing a controller... (maximize pessimistic)"
    res_mat = run_unbounded_imdp_synthesis(pimdp_filename)

    @info "Synthesizing a controller... (maximize optimistic)"
    res_mat_opt = run_imdp_synthesis(pimdp_filename, -1, mode2="optimistic", save_mats=false)

    if add_opt
        for j in 1:length(res_mat[:,1])
            if res_mat[j, 3] == res_mat_opt[j,3] && res_mat_opt[j,4] > res_mat[j,4] 
                res_mat[j, 2] = res_mat_opt[j,2]
                res_mat[j, 4] = res_mat_opt[j,4] 
            end
        end
    end

    # Copy region data
    sys_dir = res_dirs[1]
    cp("$sys_dir/regions.bson", "$dst_dir/regions.bson", force=true)
    plot_gamma_value(dst_dir, res_mat, res_mat_opt, num_dfa_states=length(dfa.states))
    save_legacy_mats(res_mat, dst_dir, -1)
    
    # plot_results_from_file(dst_dir)
    plot_synthesis_results(dst_dir, res_mat, imdp, dfa, pimdp)

    if !isnothing(modes)
        @info "Running a simulation using the resulting controller..."
        trajectories = []
        mt = MersenneTwister(111)

        if !isnothing(sim_points)
            num_sims=length(sim_points)
        end

        for i=1:num_sims
            if !isnothing(sim_points)
                x0 = sim_points[i]
            else    
                x0 = 4*rand(mt,2,1)[:]-[2.; 2.]
            end
            res = BSON.load(regions_file)
            noise_dist = TruncatedNormal(0, 0.01, -0.01, 0.01)
            simulate_system(x0, modes, 100, imdp, pimdp, dfa, res[:extents], res_mat[:,2], noise_dist=noise_dist)
            push!(trajectories, copy(pimdp.trajectory))
        end
        plot_synthesis_results(dst_dir, res_mat, imdp, dfa, pimdp, trajectories=trajectories, filename="sim.png")
    end

    return res_mat, dst_dir, pimdp
    # Plot Results
end

"
Creates the experiment directory form the parameter structure.
"
function create_experiment_directory(params)
    if !isnothing(params.data_params)
        data_tag = @sprintf("m%d-σ%1.3f-rs%d", params.data_params.data_num, params.data_params.noise_sigma, params.random_seed)
    else
        data_tag = "known-system"
    end
    for dim_key in keys(params.discretization_step)
        data_tag = @sprintf("%s-%0.3f", data_tag, params.discretization_step[dim_key])
    end
    exp_dir = @sprintf("%s/modes/%s/%s", params.experiment_directory, params.system_params.mode_tag, data_tag)
    !isdir(exp_dir) && mkpath(exp_dir)
    return exp_dir
end

function initialize_log(params; logging=Logging.Info)
    exp_dir = create_experiment_directory(params)
    # Setup the logfile with the desired level.
    glogger = SimpleLogger(stdout, logging)
    logfile = open("$exp_dir/log.txt", "w+")
    text_logger = SimpleLogger(logfile, logging)
    demux_logger = TeeLogger(glogger, text_logger)
    global_logger(demux_logger)
    @info "Experiment directory: " exp_dir
    return logfile
end

function generate_estimates(params, x_train, y_train; reuse_gp_flag=false)
    exp_dir = create_experiment_directory(params)
    # Generate a set of GPRs if none are provided.
    gps_dir = @sprintf("%s/gps", params.experiment_directory)
    !isdir(gps_dir) && mkpath(gps_dir)
    gps_filename = @sprintf("%s/%s-m%d-σ%1.3f-rs%d-gps.bin", gps_dir, params.system_params.mode_tag, params.data_params.data_num, params.data_params.noise_sigma, params.random_seed)
    data_deps = params.system_params.dependency_dims
    if isfile(gps_filename) && reuse_gp_flag
        @info "Resuing gp: " gps_filename
        open(gps_filename) do f
            gp_set = deserialize(f)
        end
    else 
        # Train GPs on this system
        gp_set = Dict()
        # TODO: This is not hyperparameterized.
        ls = 0.65
        for (i, out_dim) in enumerate(keys(data_deps)) 
            # Handle data dependency here
            x_train_sub = x_train[:, findall(.>(0), data_deps[out_dim])[:]]
            m_prior = MeanZero()
            k_prior = SE(ls, 0.)
            lnoise = log(sqrt(1+2/params.data_params.data_num)) # Generalize to handle any bound
            opt_idx = StatsBase.sample(1:length(y_train[:,1]), params.data_params.data_num_optimize, replace = false)
            gp_pre = GP(x_train_sub[opt_idx, :]', y_train[opt_idx,i], m_prior, k_prior, lnoise) 
            optimize!(gp_pre)
            gp = GP(x_train_sub', y_train[:,i], m_prior, k_prior, lnoise)
            gp_set["x$i"] = deepcopy(gp)
        end
    end

    gp_info_dict = Dict()
    for dim_key in keys(data_deps) 
        gp_info_dict[dim_key] = create_gp_info(params, gp_set, dim_key) 
    end

    # TODO: Plotting should be exposed to get accurate timing
    # create_plots ? plot_gp_fields(exp_dir, dyn_fn) : nothing

    return gp_set, gp_info_dict
end

function save_gp_info(params, gp_set, gp_info_dict)
    exp_dir = create_experiment_directory(params)
    gps_dir = @sprintf("%s/gps", params.experiment_directory)
    !isdir(gps_dir) && mkpath(gps_dir)
    gps_filename = @sprintf("%s/%s-m%d-σ%1.3f-rs%d-gps.bin", gps_dir, params.system_params.mode_tag, params.data_params.data_num, params.data_params.noise_sigma, params.random_seed)
    # Save the GPs for further analysis. 
    @info "Saving GP regressions to experiment directory..."
    open(gps_filename, "w") do f
        serialize(f, gp_set)
    end
    local_gp_file = "$exp_dir/gaussian_processes_set.bin"
    open(local_gp_file, "w") do f
        serialize(f, gp_set)
    end
    # TODO: Save the entire info dict
end

function generate_region_images(params, gp_info_dict; reuse_regions_flag=false)

    exp_dir = create_experiment_directory(params)
    region_filename = "$exp_dir/regions.bson"

    if reuse_regions_flag && isfile(region_filename)
        region_data = BSON.load(region_filename)
    else
        region_dict, region_pairs = create_region_data(params.domain, params.discretization_step)
        num_regions = length(keys(region_dict))
        region_post_dict = Dict()

        # Minimum and maximum extents
        domain = params.domain
        min_ex = [domain[dim_key][1] for dim_key in keys(domain)]
        max_ex = [domain[dim_key][2] for dim_key in keys(domain)]

        # σ_U_ubs = []
        # for dim_key in keys(domain)
        #     _, σ_U_lb, σ_U_ub = compute_σ_ub_bounds_auto(gp_info_dict[dim_key].gp, min_ex, max_ex) 
        #     push!(σ_U_ubs, σ_U_ub)
        # end
        # The last state!!
        # σ_bound_dict[num_regions] = σ_U_ubs 

        #Try to parallelize this loop.

        # Store the predicted mean and covariance for each sampled point 
        dim_keys = keys(domain)
        for i=1:num_regions-1
            @info "Bounding region $i/$num_regions"
            extent = region_dict[i]
            lb = [extent[dim_key][1] for dim_key in dim_keys]
            ub = [extent[dim_key][2] for dim_key in dim_keys]
            region_post_dict[i] = bound_extent(extent, lb, ub, gp_info_dict, dim_keys, params.system_params.dependency_dims; known_part_flag=!isnothing(params.system_params.known_dynamics_fcn))
        end  # End threaded forloop

        region_data = Dict()
        region_data[:pairs] = region_pairs
        region_data[:extents] = region_dict
        region_data[:posts] = region_post_dict
    end

    return region_data
end

function save_region_data(params, region_data)
    @info "Saving the region data to the experiment directory..."
    exp_dir = create_experiment_directory(params)
    region_filename = "$exp_dir/regions.bson"

    bson(region_filename, region_data)
end

"Generates transition bounds with the given parameters."
function generate_transition_bounds(params, gp_info_dict, region_data; 
                                    reuse_mats_flag=false) 

    exp_dir = create_experiment_directory(params)
    region_pairs = region_data[:pairs]
    region_dict = region_data[:extents]
    region_post_dict = region_data[:posts]

    tmat_filename = "$exp_dir/transition_mats.mat"
    if !reuse_mats_flag || !isfile(tmat_filename)
        @info "Calculating transition probability bounds between regions..."
        results_df = process(region_pairs, region_dict, region_post_dict, params.data_params.epsilon, gp_info_dict, params)

        # Save the matrices for use in MATLAB verification tool
        num_states = length(keys(region_dict))
        minPr_mat = zeros((num_states, num_states))
        maxPr_mat = zeros((num_states, num_states))
        for i = 1:1:num_states
            sub_row = results_df[results_df.Set1 .== i, :]
            if i == num_states
                minPr_mat[i, i] = 1.
                maxPr_mat[i, i] = 1.
            else
                for j = 1:1:num_states
                    if j == num_states
                        subsub_row = sub_row[sub_row.Set2 .== -11, :]
                        minPr_mat[i, j] = 1. - subsub_row.MaxPr[1]
                        maxPr_mat[i, j] = 1. - subsub_row.MinPr[1]
                    else
                        subsub_row = sub_row[sub_row.Set2 .== j, :]
                        minPr_mat[i, j] = subsub_row.MinPr[1]
                        maxPr_mat[i, j] = subsub_row.MaxPr[1]
                    end
                end
            end
        end
    end

    res_mats = Dict("minPr" => minPr_mat, "maxPr" => maxPr_mat)
    return res_mats 
end

function save_transition_matrices(params, res_mats)
    exp_dir = create_experiment_directory(params)
    tmat_filename = "$exp_dir/transition_mats.mat"
    matwrite(tmat_filename, res_mats) 
end

function perform_imdp_verification(params, res_mats, horizon; simply_safe=true)
    exp_dir = create_experiment_directory(params)
    imdp = create_single_mode_imdp(res_mats["minPr"], res_mats["maxPr"])
    # TODO: Need to generalize this!!!!!!
    if simply_safe
        unsafe_states = [length(imdp.states)]
        # All of this should be inside "verify IMDP"
        pimdp = construct_PIMDP_from_IMDP(imdp, unsafe_states)
        filename = "$exp_dir/test-pimdp.txt"
        write_pimdp_to_file(pimdp, filename)
        result_mat = run_bounded_imdp_verification(filename, horizon)
        safety_result_mat = zeros(size(result_mat))
        safety_result_mat[:, 1] = result_mat[:, 1]
        safety_result_mat[:, 2] = result_mat[:, 2]
        safety_result_mat[:, 3] = 1 .- result_mat[:, 4]
        safety_result_mat[:, 4] = 1 .- result_mat[:, 3]
    end

    save_legacy_mats(safety_result_mat, exp_dir, horizon)
    return safety_result_mat
end

function save_time_info(params, time_info)
    exp_dir = create_experiment_directory(params)
    time_filename = "$exp_dir/time_info.bson"
    bson(time_filename, time_info)
end

function end_to_end_transition_bounds(params; single_mode_verification=false, reuse_regions_flag=false) 
    logfile = initialize_log(params)
    @info "Generating the training data..."
    total_runtime = 0.
    timing_info = Dict()
    gen_time = @elapsed begin
        x_train, y_train = generate_training_data(params)
    end
    total_runtime += gen_time
    timing_info["data_generation_time_s"] = gen_time
    @info "Data generation time: " gen_time
    @info "Generating the regressions..."
    gp_time = @elapsed begin
        gp_set, gp_info_dict = generate_estimates(params, x_train, y_train)
    end
    total_runtime += gp_time
    timing_info["gp_regression_time_s"] = gp_time 
    @info "GP generation time: " gp_time
    save_gp_info(params, gp_set, gp_info_dict)
    @info "Generating the region info..."
    region_time = @elapsed begin
        region_data = generate_region_images(params, gp_info_dict, reuse_regions_flag=reuse_regions_flag)
    end
    total_runtime += region_time
    timing_info["region_bound_time_s"] = region_time 
    @info "Region generation time: " region_time
    save_region_data(params, region_data)
    @info "Generating the transition bounds..."
    bound_time = @elapsed begin
        result_mats = generate_transition_bounds(params, gp_info_dict, region_data)
    end
    total_runtime += bound_time
    timing_info["transition_bound_time_s"] = bound_time 
    @info "Bound generation time: " bound_time
    save_transition_matrices(params, result_mats)
 
    safety_result_mat = nothing
    if single_mode_verification
        @info "Performing safety verification on single mode..."
        verification_time = @elapsed begin 
            # script_path = "/Users/john/Projects/Julia/BMDPVerification"
            # run_MATLAB_verification(script_path, create_experiment_directory(params), -1, params.verification_mode)
            safety_result_mat = perform_imdp_verification(params, result_mats, 1)
        end
        timing_info["single_verification_time_s"] = verification_time  
        @info "Verification time: " verification_time
    end

    @info "Total runtime: " total_runtime
    timing_info["total_runtime_s"] = total_runtime
    save_time_info(params, timing_info)
    
    flush(logfile)
    close(logfile)

    return timing_info, safety_result_mat, result_mats
end

###############################################################################################
# Local Methods
###############################################################################################

function create_gp_info(params, gp_set, dim_key)
    gp = gp_set[dim_key]

    scale_factor = params.data_params.bound_type == "rkhs-tight" ? params.data_params.noise_sigma/sqrt(1. + 2. / params.data_params.data_num) : 1.
    @info "Scale factor: $scale_factor"
    domain = params.domain
    diam_domain = 0.
    for dim_key in keys(domain)
        diam_domain += (domain[dim_key][1] - domain[dim_key][2])^2
    end
    diam_domain = sqrt(diam_domain)
    σ_inf = sqrt(gp.kernel.σ2*exp(-1/2*(diam_domain)^2/gp.kernel.ℓ2))

    # Calculating the RKHS parameter bounds
    RKHS_bound = abs(domain[dim_key][2] + params.system_params.lipschitz_bound*diam_domain)/ σ_inf
    @info "RKHS Norm Bound in $dim_key: ", RKHS_bound

    B = 1 + (params.data_params.noise_sigma)^(-2)
    γ = 0.5*params.data_params.data_num*log(B)
    @info "Info gain term in $dim_key: ", γ

    K_inv = inv(gp.cK.mat + exp(gp.logNoise.value)^2*I)
    gp_info = GPInfo(gp, γ, RKHS_bound, params.data_params.bound_type, scale_factor, K_inv)

    return gp_info
end

# TODO: Fully implement this
function bound_extent(extent, lb, ub, gp_info_dict, dim_keys, data_deps; known_part_flag=false)

    lb = [extent[dim_key][1] for dim_key in dim_keys]
    ub = [extent[dim_key][2] for dim_key in dim_keys] 
    post_extent = Dict()
    σ_bounds = []

    for dim_key in dim_keys
        lbf = lb[findall(.>(0), data_deps[dim_key][:])]
        ubf = ub[findall(.>(0), data_deps[dim_key][:])]
        x_lb, μ_L_lb, μ_L_ub = compute_μ_bounds_bnb(gp_info_dict[dim_key].gp, lbf, ubf) 
        x_ub, μ_U_lb, μ_U_ub = compute_μ_bounds_bnb(gp_info_dict[dim_key].gp, lbf, ubf, max_flag=true)
        _, σ_U_lb, σ_U_ub = compute_σ_ub_bounds(gp_info_dict[dim_key].gp, gp_info_dict[dim_key].Kinv, lbf, ubf) 
        @assert μ_L_lb <= μ_U_ub
        if !(σ_U_lb <= σ_U_ub)
            @error "Sigma LB is greater than sigma UB! ", σ_U_lb, σ_U_ub
            throw("aborting")
        end 
        if known_part_flag
            # Assume an identity form for now. Requires its own bounding process.
            post_extent[dim_key] = [extent[dim_key][1] + μ_L_lb, extent[dim_key][2] + μ_U_ub] 
        else
            post_extent[dim_key] = [μ_L_lb, μ_U_ub]
        end
        push!(σ_bounds, σ_U_ub)
    end
    return post_extent, σ_bounds
end

function fetch_results(N, r)
    results_df = DataFrame(Set1 = Int[], Set2 = Int[], MinPr = Float64[], MaxPr = Float64[], MinPrPoint = Array[], MaxPrPoint = Array[])
    for j in 1:N
        total_results = fetch(r[j])
        push!(results_df, total_results)
    end
    return results_df 
end

function process(region_pairs, region_dict, region_post_dict, epsilon, gp_dict, params)
    N = length(region_pairs)
    r = Array{Array, 1}(undef, N)

    for i=1:N
        r[i] = region_pair_transitions(region_pairs[i], region_dict, region_post_dict, epsilon, gp_dict, params)
    end

    results_df = fetch_results(N, r)
    return results_df
end

function create_region_data(space, grid_sizes)
    region_dict = discretize_set(space, grid_sizes)

    # THIS IS NOT CORRECT for a 3d model. 
    region_dict[-11] = space # Region with index -11 corresponds to the safe set as a whole. 
    # region_dict[-11]["x3"] = [-1e4, 1e4]

    # Construct a list of all region pairs
    starting_regions = collect(1:1:length(keys(region_dict))-1)  
    region_pairs = []
    for starting_region in starting_regions
        for target_region_id in keys(region_dict) 
            push!(region_pairs, (starting_region, target_region_id))
        end
    end   
    return region_dict, region_pairs
end

function create_region_data_new(space, grid_size)
    region_dict, extent_dict = discretize_set_lazy(space, grid_size)
    x = space["x1"]
    y = space["x2"]
    # TODO: How to deal with this negative index in a better way
    region_dict[-11] = VPolygon([[x[1], y[1]], [x[1], y[2]], 
                                    [x[2], y[2]], [x[2], y[1]]])
    extent_dict[-11] = space

    # Construct a list of all region pairs
    starting_regions = collect(1:1:length(keys(region_dict))-1)  
    region_pairs = []
    for starting_region in starting_regions
        for target_region_id in keys(region_dict) 
            push!(region_pairs, (starting_region, target_region_id))
        end
    end   
    return region_dict, region_pairs, extent_dict
end

end
