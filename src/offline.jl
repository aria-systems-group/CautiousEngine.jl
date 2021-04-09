# End-to-end methods for offline verification and synthesis


function perform_synthesis_from_result_dirs(res_dirs, exp_dir, system_tag, spec_file, label_fcn; 
                                            plot_graphs=true, rerun_flag=false, modes=nothing, add_opt=false,
                                            region_specific="", num_sims=10, sim_points=nothing)

    @info "Creating the IMDP..."
    imdp = create_imdp_from_result_dirs(res_dirs, "$exp_dir/$system_tag")    
    regions_file = @sprintf("%s/regions%s.bson", res_dirs[1], region_specific)
    create_imdp_labels(label_fcn, imdp, regions_file)
    # imdp.labels[length(keys(imdp.labels))] = "!a∧b"
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
    res = BSON.load(regions_file)
    simulate_system(x0, modes, 100, imdp, pimdp, dfa, res[:extents], res_mat[:,2])
    push!(trajectories, copy(pimdp.trajectory))
    end
    plot_synthesis_results(dst_dir, res_mat, imdp, dfa, pimdp, trajectories=trajectories, filename="sim.png")
    end

    return res_mat, dst_dir, pimdp
    # Plot Results
end