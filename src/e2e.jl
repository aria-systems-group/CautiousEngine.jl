#== What I want from an End to end function:
3. Save the experimental parameters (as a serializable object?)
4. For each component, save the time and a sparse rep. of the result
5. An easy to plot visualization, separate from these functions
==#

"""
Performs k-step refinement over uncertain states according to safety verification.
"""
function safety_based_refinement(experiment_params::ExperimentParameters, results_path::String, verification_mat, refinement_steps::Int; 
                                 safety_threshold=0.95, horizon=1, reuse_regions_flag=false, reuse_transition_mats_flag=false)

    working_dir = results_path
    refinement_result_dirs = []

    for i=1:refinement_steps
        results_dir = "$results_path/refinement$i"
        states_to_update = findall(x->x<safety_threshold, verification_mat[1:end,3]) ∩ findall(x->x>safety_threshold, verification_mat[1:end,4])  
        imdp = refine_imdp(states_to_update, experiment_params, working_dir=working_dir, results_dir=results_dir, reuse_regions_flag=reuse_regions_flag, reuse_transition_mats_flag=reuse_transition_mats_flag)
        verification_mat = Globally(imdp, "safe", horizon, "$results_dir/imdp-r$i.txt")
        save_legacy_mats(verification_mat, results_dir, horizon)
        working_dir = results_dir 
        push!(refinement_result_dirs, results_dir)
    end

    return refinement_result_dirs
end

function safety_based_synthesis(experiment_params_array, system_paths::Array{String,1}, results_path::String, 
    system_tag::String, label_fcn, refinement_steps::Int; 
    minimum_threshold=0.80, horizon=-1, reuse_regions_flag=false, reuse_transition_mats_flag=false) 

    exp_dir = results_path

    synth_dir = "$results_path/safety"
    synthesis_mat, results_path, _ = perform_synthesis_from_result_dirs_safety(system_paths, synth_dir, system_tag, label_fcn); 
    save_legacy_mats(synthesis_mat, synth_dir, horizon)

    refinement_result_dirs = []

    for i=1:refinement_steps
        results_dir = "$synth_dir/refinement$i"
        states_to_update = findall(x->x<minimum_threshold, synthesis_mat[1:end,3]) ∩ findall(x->x>minimum_threshold, synthesis_mat[1:end,4])  

        new_imdp_dirs = Array{String,1}()
        # For each IMDP in the system
        for (j, component_dir) in enumerate(system_paths)
            imdp_results_dir = "$results_dir/mode$j"
            refine_imdp(states_to_update, experiment_params_array[j], working_dir=component_dir, results_dir=imdp_results_dir, reuse_regions_flag=reuse_regions_flag, reuse_transition_mats_flag=reuse_transition_mats_flag)
            push!(new_imdp_dirs, imdp_results_dir)
        end

        cp("$results_dir/mode1/regions.bson", "$results_dir/regions.bson", force=true)

        synthesis_mat, _, _ = perform_synthesis_from_result_dirs_safety(new_imdp_dirs, results_dir, system_tag, label_fcn; horizon=horizon); 

        save_legacy_mats(synthesis_mat, results_dir, horizon)

        system_paths = new_imdp_dirs
        push!(refinement_result_dirs, results_dir)
    end

    return refinement_result_dirs
end

"""
Performs k-step refinement over uncertain states according to minimum spec satisfaction from synthesis.
"""
function spec_based_refinement(experiment_params_array, system_paths::Array{String,1}, results_path::String, 
                               system_tag::String, label_fcn, refinement_steps::Int; 
                               minimum_threshold=0.80, horizon=-1, reuse_regions_flag=false, reuse_transition_mats_flag=false)

    
    refinement_result_dirs = []
    synth_dir = "$results_path/$system_tag/reach-specification"
    !isdir(synth_dir) && mkpath(synth_dir)

    res_mat, res_dir, pimdp = perform_synthesis_from_result_dirs(system_paths, synth_dir, "", experiment_params_array[1].specification_file, label_fcn; add_opt=true, horizon=horizon)

    save_legacy_mats(res_mat, synth_dir, horizon)

    for i=1:refinement_steps
        results_dir = "$synth_dir/refinement$i"
        states_to_update = unique(Int.(pimdp_col_to_imdp_state.(findall(x->x<minimum_threshold, res_mat[1:end-1,3]) ∩ findall(x->x>minimum_threshold, res_mat[1:end-1,4]), 2))) # ! TODO: FIX THIS

        new_imdp_dirs = Array{String,1}()
        # For each IMDP in the system
        for (j, component_dir) in enumerate(system_paths)
            imdp_results_dir = "$results_dir/mode$j"
            refine_imdp(states_to_update, experiment_params_array[j], working_dir=component_dir, results_dir=imdp_results_dir, reuse_regions_flag=reuse_regions_flag, reuse_transition_mats_flag=reuse_transition_mats_flag)
            push!(new_imdp_dirs, imdp_results_dir)
        end
       
        res_mat, res_dir, pimdp = perform_synthesis_from_result_dirs(new_imdp_dirs, results_dir, "", experiment_params_array[1].specification_file, label_fcn; add_opt=true, horizon=horizon)
        cp("$results_dir/switched-system/regions.bson", "$results_dir/regions.bson", force=true)
        system_paths = new_imdp_dirs
        push!(refinement_result_dirs, results_dir) 
    end

    return refinement_result_dirs
end

function end_to_end_transition_bounds(params; reuse_gps_flag=false, reuse_regions_flag=false) 
    
	@info "Initializing experiment and saving the parameters..."
	exp_dir = create_experiment_directory(params)
	params.experiment_directory = exp_dir
    logfile = initialize_log(exp_dir)
	save_experiment_params(params)
	total_runtime = 0.
    timing_info = Dict()

	regions_filename = "$exp_dir/regions.bson"
	gps_filename = generate_gp_filename(params)

	if reuse_regions_flag && (isfile(regions_filename) && isfile(gps_filename))
		gp_set, gp_info_dict = load_gp_info(gps_filename) 
		timing_info["data_generation_time_s"] = -1.
		timing_info["gp_regression_time_s"] = -1. 
		region_data = load_region_data("$exp_dir/regions.bson") 
		timing_info["data_generation_time_s"] = -1.
		timing_info["gp_regression_time_s"] = -1. 
		timing_info["region_bound_time_s"] = -1. 
		@goto reuse_regions
	elseif reuse_gps_flag && isfile(gps_filename)
		gp_set, gp_info_dict = load_gp_info(gps_filename) 
		timing_info["data_generation_time_s"] = -1.
		timing_info["gp_regression_time_s"] = -1. 
		@goto reuse_gps
	end

    @info "Generating the training data..."
    
    gen_time = @elapsed begin
        x_train, y_train = generate_training_data(params)
    end
    total_runtime += gen_time
    timing_info["data_generation_time_s"] = gen_time
    @info "Data generation time: " gen_time
	save_training_data(x_train, y_train, params.experiment_directory)	

    @info "Generating the GP regressions..."
    gp_time = @elapsed begin
        gp_set, gp_info_dict = generate_estimates(params, x_train, y_train)
    end
    total_runtime += gp_time
    timing_info["gp_regression_time_s"] = gp_time 
    @info "GP generation time: " gp_time
    save_gp_info(params, gp_set, gp_info_dict)

@label reuse_gps
    @info "Generating the region info..."
    region_time = @elapsed begin
        region_data = generate_region_images(params, gp_info_dict)
    end
    total_runtime += region_time
    timing_info["region_bound_time_s"] = region_time 
    @info "Region generation time: " region_time
    save_region_data(exp_dir, region_data)

@label reuse_regions
    @info "Generating the transition bounds..."
    bound_time = @elapsed begin
        result_mats = generate_transition_bounds(params, gp_info_dict, region_data)
    end
    total_runtime += bound_time
    timing_info["transition_bound_time_s"] = bound_time 
    @info "Bound generation time: " bound_time
    save_transition_matrices(params, result_mats)

    @info "Total runtime: " total_runtime
    timing_info["total_runtime_s"] = total_runtime
    save_time_info(exp_dir, timing_info)
    
    flush(logfile)
    close(logfile)

    return timing_info, result_mats
end

function end_to_end_transition_bounds_local_gps(params; single_mode_verification=false, reuse_regions_flag=false) 
    exp_dir = create_experiment_directory(params)
	params.experiment_directory = exp_dir
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

function perform_synthesis_from_result_dirs_safety(res_dirs, exp_dir, system_tag, label_fcn; 
	plot_graphs=true, rerun_flag=false, modes=nothing, add_opt=true,
	region_specific="", num_sims=10, sim_points=nothing, horizon=-1)

    @info "Creating the IMDP..."
	imdp = create_imdp_from_result_dirs(res_dirs, "$exp_dir")    
    @info("$exp_dir")
	regions_file = @sprintf("%s/regions%s.bson", res_dirs[1], region_specific)
	create_imdp_labels(label_fcn, imdp, regions_file)
	plot_graphs ? create_dot_graph(imdp, "$exp_dir/imdp.gv") : nothing

	# res_mat = run_pimdp_synthesis(pimdp, pimdp_filename, add_opt=add_opt)
    imdp_filename = "$exp_dir/imdp.txt"
    res_mat = CautiousEngine.Globally(imdp, "safe", horizon, imdp_filename, synthesis_flag=true)

    if add_opt
        @info "Synthesizing a controller... (maximize optimistic)"
        res_mat_opt = run_imdp_synthesis(imdp_filename, horizon, mode2="optimistic", save_mats=false)
        for j in 1:length(res_mat[:,1])
            if res_mat[j, 3] == res_mat_opt[j,3] && res_mat_opt[j,4] > res_mat[j,4] 
                res_mat[j, 2] = res_mat_opt[j,2]
                res_mat[j, 4] = res_mat_opt[j,4] 
            end
        end
    end

    dst_dir = "$exp_dir/safety"
	isdir(dst_dir) ? nothing : mkpath(dst_dir)
	# Copy region data
	sys_dir = res_dirs[1]
	cp("$sys_dir/regions.bson", "$dst_dir/regions.bson", force=true)

	# TODO Plot gamma values somewhere else, i.e. not here
	# if add_opt
	#     plot_gamma_value(dst_dir, res_mat, res_mat_opt, num_dfa_states=length(dfa.states))
	# end
	save_legacy_mats(res_mat, dst_dir, horizon)

	# plot_results_from_file(dst_dir)
	# plot_synthesis_results(dst_dir, res_mat, imdp, dfa, pimdp)
    return res_mat, dst_dir, -1 
end

function perform_synthesis_from_result_dirs(res_dirs, exp_dir, system_tag, spec_file, label_fcn; 
	plot_graphs=false, rerun_flag=false, modes=nothing, add_opt=false,
	region_specific="", num_sims=10, sim_points=nothing, horizon=-1)

	@info "Creating the IMDP..."
	imdp = create_imdp_from_result_dirs(res_dirs, "$exp_dir/$system_tag")    
	regions_file = @sprintf("%s/regions%s.bson", res_dirs[1], region_specific)
	create_imdp_labels(label_fcn, imdp, regions_file)
	plot_graphs ? create_dot_graph(imdp, "$exp_dir/$system_tag/imdp.gv") : nothing

	@info "Creating the DFA..."
	spec_tag = basename(spec_file)[1:end-3]
	include("$spec_file")
	dfa = DFA(dfa_states, dfa_props, dfa_transitions, dfa_accepting_state, dfa_sink_state, dfa_initial_state)
	dst_dir = "$exp_dir/$system_tag/$spec_tag"
	isdir(dst_dir) ? nothing : mkpath(dst_dir)
	plot_graphs ? create_dot_graph(dfa, "$dst_dir/specification.gv") : nothing

	@info "Creating pimdp..."
	pimdp = construct_DFA_IMDP_product(dfa, imdp)
	plot_graphs ? create_dot_graph(pimdp, "$dst_dir/pimdp.gv") : nothing

	pimdp_filename = "$dst_dir/pimdp.txt"
	res_mat = run_pimdp_synthesis(pimdp, pimdp_filename, add_opt=add_opt, horizon=horizon)

	# Copy region data
	sys_dir = res_dirs[1]
	cp("$sys_dir/regions.bson", "$dst_dir/regions.bson", force=true)

	# TODO Plot gamma values somewhere else, i.e. not here
	# if add_opt
	#     plot_gamma_value(dst_dir, res_mat, res_mat_opt, num_dfa_states=length(dfa.states))
	# end
	save_legacy_mats(res_mat, dst_dir, horizon)

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
			extents, _ = deserialize_region_data(regions_file)
			simulate_system(x0, modes, 100, imdp, pimdp, dfa, extents, res_mat[:,2])
			push!(trajectories, copy(pimdp.trajectory))
		end
		plot_synthesis_results(dst_dir, res_mat, imdp, dfa, pimdp, trajectories=trajectories, filename="sim.png")
	end

	return res_mat, dst_dir, pimdp
end