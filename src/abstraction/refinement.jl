"""
Refines the specificed IMDP states using a uniform discretization of each state.
# Arguments
- `states_to_refine` - Vector of state indeces to be refined
- `params::ExperimentParameters` - Experiment parameters structure
"""
function refine_imdp(states_to_refine, params::ExperimentParameters; 
    working_dir=nothing, results_dir=nothing, 
    reuse_regions_flag=false, reuse_transition_mats_flag=false,
    label_fcn=nothing)

    experiment_dir = params.experiment_directory
    if isnothing(working_dir)
        working_dir = experiment_dir
    end

    if isnothing(results_dir)
        results_dir = "$exp_dir/refined"
    end
    !isdir(results_dir) && mkdir(results_dir)

    tmat_filename = "$results_dir/transition_mats.bson"
    # ! Transition matrices being saved incorrectly. Not able to reload.
    if reuse_transition_mats_flag && isfile(tmat_filename) && false
        transition_mats = load_transition_matrices(tmat_filename)
        @goto reuse_transition_mats
    end

    # Load the GPs
    gp_filename = generate_gp_filename(params)
    _, gp_info = load_gp_info(gp_filename)

    region_filename = "$results_dir/regions.bson"
    if reuse_regions_flag && isfile(region_filename)
        region_data = load_region_data(region_filename)
        @goto reuse_regions
    end

    region_data = load_region_data("$working_dir/regions.bson") 
    extents = region_data[:extents]
    post_data = region_data[:posts]

    # Generate new discretization
    new_region_ids = []
    for state_to_refine in states_to_refine
        extents, new_ids = uniform_refinement(state_to_refine, extents)
        new_region_ids = new_region_ids ∪ new_ids 
    end

    # Generate new posterior bounds
    Threads.@threads for new_id in new_region_ids
        # TODO: Add a check here for "small enough" covariance bounds, and just update mean if is
        # TODO: Add way to do the fast extent bounding
        post_data[new_id] = bound_extent(extents[new_id], gp_info, params.system_params.dependency_dims)
    end

    # Fix the extent ids to correspond to order
    new_extent_dict = Dict()
    new_post_dict = Dict()
    sorted_keys = sort(collect(keys(extents)))
    for (i,j) in zip(1:extents.count-1, sorted_keys[1:end-1])
        new_extent_dict[i] = extents[j]
        new_post_dict[i] = post_data[j] 
    end
    new_extent_dict[new_extent_dict.count+1] = extents[sorted_keys[end]]

    # Overwrite region data dictionary 
    region_data[:extents] = new_extent_dict
    region_data[:posts] = new_post_dict
    save_region_data(new_extent_dict, new_post_dict, region_filename)

    @label reuse_regions
    transition_mats = generate_transition_bounds(params, gp_info, region_data)
    save_transition_matrices(transition_mats, tmat_filename)

    @label reuse_transition_mats
    # Rebuild the IMDP through rebuilding the matrices
    new_imdp = create_simple_imdp(transition_mats["minPr"], transition_mats["maxPr"])

    if !isnothing(label_fcn)
        create_imdp_labels(label_fcn, new_imdp, new_extent_dict)
    end

    return new_imdp
end