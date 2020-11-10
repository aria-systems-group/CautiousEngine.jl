function create_switched_system(system_dirs, exp_dir)

    # Create a new system dir with the tags
    new_dirname = "$exp_dir/switched-system" 
    !isdir(new_dirname) && mkpath(new_dirname)

    # Copy the region info from the first system
    first_dir = system_dirs[1]
    cp("$first_dir/regions.bson", "$new_dirname/regions.bson", force=true)

    # Create new transition matrices consistent with imdp solver

    # Load all of the transition matrices into memory
    # Rather, take all of the imdps
    # Given a bunch of transition matrices, create the biggun

    num_modes = length(system_dirs)
    minPr_mats = []
    maxPr_mats = []

    for system_dir in system_dirs
        res = matread("$system_dir/transition_mats.mat")
        push!(minPr_mats, res["minPr"])
        push!(maxPr_mats, res["maxPr"])
    end
    mat_size = size(minPr_mats[1])[1]
    # multimode_transitions = hcat([[0.0..0.0 for i=1:mat_size*mat_size] for j=1:num_modes]...) 
    mm_minPr = hcat([[0.0 for i=1:mat_size*num_modes] for j=1:mat_size]...) 
    mm_maxPr = hcat([[0.0 for i=1:mat_size*num_modes] for j=1:mat_size]...) 

    for i=1:num_modes
        # multimode_transitions[i:num_modes:end] = transition_mats[i]
        mm_minPr[i:num_modes:end, :] = minPr_mats[i]
        mm_maxPr[i:num_modes:end, :] = maxPr_mats[i]
    end

    # minPr = minimum.(multimode_transitions)
    # maxPr = maximum.(multimode_transitions)

    # Save and return directory
    matwrite("$new_dirname/transition_mats.mat", Dict("minPr" => mm_minPr, "maxPr" => mm_maxPr))

    return new_dirname, mm_minPr, mm_maxPr 
end