"""
Creates a switched system from a multitude of single-mode systems. 

# Arguments
- `system_dirs::Array{String,1}`: directories containing the results produced by single mode systems.
- `exp_dir::String`: location to create the switched-system directory.
"""
function create_switched_system(system_dirs::Array{String,1}, exp_dir::String)

    # Create a new system dir with the tags
    new_dirname = "$exp_dir/switched-system" 
    !isdir(new_dirname) && mkpath(new_dirname)

    # Copy the region info from the first system
    first_dir = system_dirs[1]
    cp("$first_dir/regions.bson", "$new_dirname/regions.bson", force=true)

    # Load all of the transition matrices into memory
    minPr_mats = []
    maxPr_mats = []

    for system_dir in system_dirs
        res = matread("$system_dir/transition_mats.mat")
        push!(minPr_mats, res["minPr"])
        push!(maxPr_mats, res["maxPr"])
    end

    mm_minPr, mm_maxPr = create_switched_system_matrix(minPr_mats, maxPr_mats)

    # Save and return directory
    matwrite("$new_dirname/transition_mats.mat", Dict("minPr" => mm_minPr, "maxPr" => mm_maxPr))

    return new_dirname, mm_minPr, mm_maxPr 
end

"""
Creates a switched system matrix from a list of min and max probabilities of transition.
"""
function create_switched_system_matrix(minPr_mats, maxPr_mats)
    num_modes = length(minPr_mats)
    mat_size = size(minPr_mats[1], 1)
    mm_minPr = hcat([[0.0 for i=1:mat_size*num_modes] for j=1:mat_size]...) 
    mm_maxPr = hcat([[0.0 for i=1:mat_size*num_modes] for j=1:mat_size]...) 

    for i=1:num_modes
        mm_minPr[i:num_modes:end, :] = minPr_mats[i]
        mm_maxPr[i:num_modes:end, :] = maxPr_mats[i]
    end

    return mm_minPr, mm_maxPr
end