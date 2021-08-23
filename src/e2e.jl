#== What I want from an End to end function:
3. Save the experimental parameters (as a serializable object?)
4. For each component, save the time and a sparse rep. of the result
5. An easy to plot visualization, separate from these functions
==#

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

function perform_synthesis_from_result_dirs(res_dirs, exp_dir, system_tag, spec_file, label_fcn; 
	plot_graphs=true, rerun_flag=false, modes=nothing, add_opt=false,
	region_specific="", num_sims=10, sim_points=nothing)

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
	res_mat = run_pimdp_synthesis(pimdp, pimdp_filename, add_opt=add_opt)

	# Copy region data
	sys_dir = res_dirs[1]
	cp("$sys_dir/regions.bson", "$dst_dir/regions.bson", force=true)

	# TODO Plot gamma values somewhere else, i.e. not here
	# if add_opt
	#     plot_gamma_value(dst_dir, res_mat, res_mat_opt, num_dfa_states=length(dfa.states))
	# end
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
			extents, _ = deserialize_region_data(regions_file)
			simulate_system(x0, modes, 100, imdp, pimdp, dfa, extents, res_mat[:,2])
			push!(trajectories, copy(pimdp.trajectory))
		end
		plot_synthesis_results(dst_dir, res_mat, imdp, dfa, pimdp, trajectories=trajectories, filename="sim.png")
	end

	return res_mat, dst_dir, pimdp
end