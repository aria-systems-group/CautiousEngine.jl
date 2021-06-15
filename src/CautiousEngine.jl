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
using Plots

# Lazy sets, the hot new way to work with polytopes
using LazySets

# My own packages
include("GPBounding.jl")
using .GPBounding

include("transitions.jl")
include("error_bounds.jl")
include("discretization.jl")
include("data.jl")
include("verification.jl")
include("imdp.jl")
include("dfa.jl")
include("pimdp.jl")
include("multimodal.jl")
include("visualize.jl")
include("intersection.jl")
include("run_simulation.jl")
include("offline.jl")
include("online_control.jl")
include("gps.jl")

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
    local_radius::Float64
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

export SystemParameters, DataParameters, ExperimentParameters

export general_label_fcn

export perform_synthesis_from_result_dirs, generate_transition_bounds, generate_linear_truth, verify_experiment, generate_result_dir_name

export plot_results_from_file, plot_gp_field_slice

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

# TODO: Fix this function
function generate_linear_truth(params, system_matrix; single_mode_verification=false)
    logfile = initialize_log(params)
    timing_info = Dict()
    total_runtime = 0.
    # TODO: Move this
    @info "Determining the post-images of the regions under the linear map."
    region_time = @elapsed begin 
        region_dict, region_pairs, extents_dict = create_region_data_new(params.domain, params.discretization_step)

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

        # TODO: Move this constructor to its own fcn
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

    verification_result_mat = nothing
    exp_dir = create_experiment_directory(params)
    if single_mode_verification
        @info "Performing safety verification on single mode for 1 step..."
        verification_time = @elapsed begin 
            imdp = create_simple_imdp(result_mats["minPr"], result_mats["maxPr"])
            horizon=1
            verification_result_mat = Globally(imdp, "safe", horizon, "$exp_dir/imdp.txt")
            save_legacy_mats(verification_result_mat, exp_dir, horizon)
        end
        timing_info["single_verification_time_s"] = verification_time  
        @info "Verification time: " verification_time
    end

    @info "Total runtime: " total_runtime
    timing_info["total_runtime_s"] = total_runtime
    save_time_info(params, timing_info)
    
    flush(logfile)
    close(logfile)

    return timing_info, verification_result_mat, result_mats
end

# TODO: This will eventually be the end-to-end function
# function perform_synthesis_experiment(params, dynamics_fcn)
   
#     generate_transition_bounds(params, dynamics_fcn)
#     # Now do all the synthesis stuff here.
#     perform_synthesis_from_result_dirs(res_dirs, params.experiment_directory, labels, params.system_tag, params.specification_file) 

# end

# function perform_verification_from_result_dirs(res_dirs, exp_dir, system_tag, spec_file, label_fcn; 
#                                                plot_graphs=true, rerun_flag=false, modes=nothing, add_opt=false,
#                                                region_specific="", num_sims=10, sim_points=nothing)
#     @info "Creating the IMDP..."
#     imdp = create_imdp_from_result_dirs(res_dirs, "$exp_dir/$system_tag")    
#     regions_file = @sprintf("%s/regions%s.bson", res_dirs[1], region_specific)
#     create_imdp_labels(label_fcn, imdp, regions_file)
#     plot_graphs ? create_graph_from_imdp(imdp, "$exp_dir/$system_tag/imdp.gv") : nothing
#     function run_bounded_imdp_verification(imdp_file, k)
                                           
# end


"""
Run the PIMDP synthesis
"""
function run_pimdp_synthesis(pimdp::PIMDP, pimdp_filename::String; add_opt=false)
    write_pimdp_to_file(pimdp, pimdp_filename)

    @info "Synthesizing a controller... (maximize pessimistic)"
    res_mat = run_unbounded_imdp_synthesis(pimdp_filename)

    if add_opt
        @info "Synthesizing a controller... (maximize optimistic)"
        basename = pimdp_filename[1:end-3]
        opt_pimdp_filename = "$basename-optimistic.txt"
        res_mat_opt = run_imdp_synthesis(pimdp_filename, -1, mode2="optimistic", save_mats=false)
        for j in 1:length(res_mat[:,1])
            if res_mat[j, 3] == res_mat_opt[j,3] && res_mat_opt[j,4] > res_mat[j,4] 
                res_mat[j, 2] = res_mat_opt[j,2]
                res_mat[j, 4] = res_mat_opt[j,4] 
            end
        end
    end

    return res_mat
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

"""
Generate region bounds under a single GP regression.
"""
function generate_region_images(params, gp_info_dict::Dict; reuse_regions_flag=false)

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
       
        # Store the predicted mean and covariance for each sampled point 
        dim_keys = keys(domain)
        if Threads.nthreads() > 1
            Threads.@threads for i=1:num_regions-1
                # @info "Bounding region $i/$num_regions"
                extent = region_dict[i]
                lb = [extent[dim_key][1] for dim_key in dim_keys]
                ub = [extent[dim_key][2] for dim_key in dim_keys]
                region_post_dict[i] = bound_extent(extent, lb, ub, gp_info_dict, dim_keys, params.system_params.dependency_dims; known_part_flag=!isnothing(params.system_params.known_dynamics_fcn))
            end
        else
            for i=1:num_regions-1
                @info "Bounding region $i/$num_regions"
                extent = region_dict[i]
                lb = [extent[dim_key][1] for dim_key in dim_keys]
                ub = [extent[dim_key][2] for dim_key in dim_keys]
                region_post_dict[i] = bound_extent(extent, lb, ub, gp_info_dict, dim_keys, params.system_params.dependency_dims; known_part_flag=!isnothing(params.system_params.known_dynamics_fcn))
            end
        end

        region_data = Dict()
        region_data[:pairs] = region_pairs
        region_data[:extents] = region_dict
        region_data[:posts] = region_post_dict
    end

    return region_data
end

"""
Generate region images using local GP regressions for each discrete state.
"""
function generate_region_images(params, x_train, y_train; reuse_regions_flag=false)

    exp_dir = create_experiment_directory(params)
    region_filename = "$exp_dir/regions.bson"

    if reuse_regions_flag && isfile(region_filename)
        region_data = BSON.load(region_filename)
    else
        region_dict, region_pairs = create_region_data(params.domain, params.discretization_step)
        num_regions = length(keys(region_dict))
        region_post_dict = Dict()
        region_gp_dict = Dict()

        # Minimum and maximum extents
        domain = params.domain
        min_ex = [domain[dim_key][1] for dim_key in keys(domain)]
        max_ex = [domain[dim_key][2] for dim_key in keys(domain)]

        # Store the predicted mean and covariance for each sampled point 
        dim_keys = keys(domain)
        Threads.@threads for i=1:num_regions
            # Get subset of data here
            if i == num_regions
                extent = region_dict[-11]
                lb = [extent[dim_key][1] for dim_key in dim_keys]
                ub = [extent[dim_key][2] for dim_key in dim_keys]
                # Assume 2D
                gp_set, gp_info_dict = generate_estimates(params, x_train, y_train, filename_appendix=i) 
                save_gp_info(params, gp_set, gp_info_dict, save_global_gps=false, filename_appendix=i)
                region_post_dict[i] = bound_extent(extent, lb, ub, gp_info_dict, dim_keys, params.system_params.dependency_dims; known_part_flag=!isnothing(params.system_params.known_dynamics_fcn))
                region_gp_dict[i] = gp_info_dict 
            else
                extent = region_dict[i]
                lb = [extent[dim_key][1] for dim_key in dim_keys]
                ub = [extent[dim_key][2] for dim_key in dim_keys]
                # Assume 2D
                center = [mean(lb), mean(ub)]
               
                kdtree = KDTree(x_train);
                
                num_neighbors = 60;
                sub_idx, _ = knn(kdtree, center, num_neighbors, true)
                
                x_sub = @view x_train[:, sub_idx]
                y_sub = @view y_train[:, sub_idx]

                gp_set, gp_info_dict = generate_estimates(params, x_sub, y_sub, filename_appendix=i) 
                save_gp_info(params, gp_set, gp_info_dict, save_global_gps=false, filename_appendix=i)
                region_post_dict[i] = bound_extent(extent, lb, ub, gp_info_dict, dim_keys, params.system_params.dependency_dims; known_part_flag=!isnothing(params.system_params.known_dynamics_fcn))
                region_gp_dict[i] = gp_info_dict
            end
        end  # End threaded forloop

        region_data = Dict()
        region_data[:pairs] = region_pairs
        region_data[:extents] = region_dict
        region_data[:posts] = region_post_dict
        region_data[:gps] = region_gp_dict
    end

    return region_data
end

function save_region_data(params, region_data)
    @info "Saving the region data to the experiment directory..."
    exp_dir = create_experiment_directory(params)
    region_filename = "$exp_dir/regions.bson"
    region_data_save = Dict()
    region_data_save[:pairs] = region_data[:pairs]
    region_data_save[:extents] = region_data[:extents]
    region_data_save[:posts] = region_data[:posts]
    bson(region_filename, region_data_save)
end

"""
Generates transition bounds with the given parameters.
"""
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

"""
Generates transition bounds with the given parameters.
"""
function generate_transition_bounds(params, region_data; 
                                    reuse_mats_flag=false) 

    exp_dir = create_experiment_directory(params)
    region_pairs = region_data[:pairs]
    region_dict = region_data[:extents]
    region_post_dict = region_data[:posts]
    region_gp_dict = region_data[:gps]

    tmat_filename = "$exp_dir/transition_mats.mat"
    if !reuse_mats_flag || !isfile(tmat_filename)
        @info "Calculating transition probability bounds between regions..."
        results_df = process_foo(region_pairs, region_dict, region_post_dict, region_gp_dict, params.data_params.epsilon, params)

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

function save_time_info(params, time_info)
    exp_dir = create_experiment_directory(params)
    time_filename = "$exp_dir/time_info.bson"
    bson(time_filename, time_info)
end

function end_to_end_transition_bounds(params; single_mode_verification=false, reuse_regions_flag=false) 
    exp_dir = create_experiment_directory(params)
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
 
    verification_result_mat = nothing
    if single_mode_verification
        @info "Performing safety verification on single mode for 1 step..."
        verification_time = @elapsed begin 
            imdp = create_simple_imdp(result_mats["minPr"], result_mats["maxPr"])
            horizon=1
            verification_result_mat = Globally(imdp, "safe", horizon, "$exp_dir/imdp.txt")
            save_legacy_mats(verification_result_mat, exp_dir, horizon)
        end
        timing_info["single_verification_time_s"] = verification_time  
        @info "Verification time: " verification_time
    end

    @info "Total runtime: " total_runtime
    timing_info["total_runtime_s"] = total_runtime
    save_time_info(params, timing_info)
    
    flush(logfile)
    close(logfile)

    return timing_info, verification_result_mat, result_mats
end

function end_to_end_transition_bounds_local_gps(params; single_mode_verification=false, reuse_regions_flag=false) 
    exp_dir = create_experiment_directory(params)
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
    @info "Generating the regressions and region data..."
    region_time = @elapsed begin
        region_data = generate_region_images(params, x_train, y_train) 
    end
    total_runtime += region_time
    timing_info["region_bound_time_s"] = region_time 
    @info "Region generation time: " region_time
    save_region_data(params, region_data)
    @info "Generating the transition bounds..."
    bound_time = @elapsed begin
        result_mats = generate_transition_bounds(params, region_data)
    end
    total_runtime += bound_time
    timing_info["transition_bound_time_s"] = bound_time 
    @info "Bound generation time: " bound_time
    save_transition_matrices(params, result_mats)
 
    safety_result_mat = nothing
    if single_mode_verification
        @info "Performing safety verification on single mode for 1 step..."
        verification_time = @elapsed begin 
            imdp = create_simple_imdp(result_mats["minPr"], result_mats["maxPr"])
            horizon=1
            verification_result_mat = Globally(imdp, "safe", horizon, "$exp_dir/imdp.txt")
            save_legacy_mats(verification_result_mat, exp_dir, horizon)
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

function process_foo(region_pairs, region_dict, region_post_dict, region_gp_dict, epsilon, params)
    N = length(region_pairs)
    r = Array{Array, 1}(undef, N)

    for i=1:N
        r[i] = region_pair_transitions(region_pairs[i], region_dict, region_post_dict, epsilon, region_gp_dict[region_pairs[i][1]], params)
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
