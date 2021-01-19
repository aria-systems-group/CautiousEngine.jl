function BoundedUntil(imdp, phi1, phi2, k, imdp_filepath)

    labels_vector = Array{String}(undef, length(imdp.states))
    for label_key in keys(imdp.labels)
        labels_vector[label_key] = imdp.labels[label_key]
    end
    
    # Get the Qyes and Qno states
    Qyes = findall(x->x==phi2, labels_vector)
    @info Qyes
    Qno = isnothing(phi1) ? [] : findall(x->!(x==phi1 || x==phi2), labels_vector)
    @info Qno

    # Write them to the file
    write_imdp_to_file_bounded(imdp, Qyes, Qno, imdp_filepath)

    # Do verification 
    if !isnothing(phi1)
        result_mat = run_bounded_imdp_verification(imdp_filepath, k)
    else
        result_mat = run_imdp_synthesis(imdp_filepath, k; ep=1e-6, mode1="maximize", mode2="optimistic", save_mats=true)
    end
        # Done-zoe
    return result_mat 
end

function Globally(imdp, phi, k, imdp_filepath)

    phi1 = nothing
    phi2 = "!$phi"

    result_mat = BoundedUntil(imdp, phi1, phi2, k, imdp_filepath)
    safety_result_mat = zeros(size(result_mat))
    safety_result_mat[:, 1] = result_mat[:, 1]
    safety_result_mat[:, 2] = result_mat[:, 2]
    safety_result_mat[:, 3] = 1 .- result_mat[:, 4]
    safety_result_mat[:, 4] = 1 .- result_mat[:, 3]
    return safety_result_mat
end

"""
Call the synthesis tool with the given parameters.
"""
function run_imdp_synthesis(imdp_file, k; ep=1e-6, mode1="maximize", mode2="pessimistic", save_mats=true)
    exe_path = "synthesis"  # Assumes that this program is on the user's path
    res = read(`$exe_path $mode1 $mode2 $k 0.000001 $imdp_file`, String)
    filter_res = replace(res, "\n"=>" ")
    res_split = split(filter_res)
    dst_dir = dirname(imdp_file)
    open("$dst_dir/verification-result-$k.txt", "w") do f
        print(f, res) 
    end
    print(res)
    res_mat = res_to_numbers(res)
    return res_mat
end

"""
Call the synthesis tool with unbounded time horizon.
"""
function run_unbounded_imdp_synthesis(imdp_file)
    res_mat = run_imdp_synthesis(imdp_file, -1)
    return res_mat
end

"""
Call the synthesis tool to perform verification on a bounded horizon.
"""
function run_bounded_imdp_verification(imdp_file, k)
    res_mat = run_imdp_synthesis(imdp_file, k, mode1="minimize", mode2="pessimistic")
    return res_mat
end

"""
Call the MATLAB verification tool.
"""
function run_MATLAB_verification(script_path, exp_dir, k, verification_mode)
    @info "Performing BMDP Verification..."
    # mat"cd($script_path)"
    # mat"addpath(genpath('./'))"
    # mat"generateVerification($exp_dir,$verification_mode,$k)"
end

"""
Call the synthesis tool to perform verification on an unbounded horizon.
"""
function run_unbounded_imdp_verification(imdp_file)
    res_mat = run_imdp_synthesis(imdp_file, -1, mode1="minimize", mode2="pessimistic")
    return res_mat
end

"""
Confirm that the verification of a data-driven system is consistent with the known system.
"""
function verify_experiment(path_to_experiment, path_to_truth)

    @info "Verifying experiment at $path_to_experiment"
    # load verification matrices
    exp_mats = MAT.matread("$path_to_experiment/globally-safe-1-steps.mat")
    # load truth matrices
    true_mats = MAT.matread("$path_to_truth/globally-safe-1-steps.mat")

    if sum(exp_mats["indVmin"] .> true_mats["indVmin"]) != 0 
        @error "Error in the following regions: " findall(x->x, exp_mats["indVmin"] .> true_mats["indVmin"]) 
    end
    @assert sum(exp_mats["indVmax"] .< true_mats["indVmax"]) == 0 

end

"""
Creates a matrix from the string output of the BMDP synthesis tool.
"""
function res_to_numbers(res_string)
    filter_res = replace(res_string, "\n"=>" ")
    res_split = split(filter_res)
    num_rows = Int(length(res_split)/4)

    res_mat = zeros(num_rows, 4)
    for i=1:num_rows
        res_mat[i, 1] = parse(Int, res_split[(i-1)*4+1])+1.
        res_mat[i, 2] = parse(Int, res_split[(i-1)*4+2])+1.
        res_mat[i, 3] = parse(Float64, res_split[(i-1)*4+3])
        res_mat[i, 4] = parse(Float64, res_split[(i-1)*4+4])
    end

    return res_mat
end

"""
Saves the output of the synthesis tool in a legacy way.
"""
function save_legacy_mats(res_mat, dst_dir, k)
    policy = res_mat[:, 2]
    indVmax = res_mat[:, 4]
    indVmin = res_mat[:, 3]

    filename = "$dst_dir/globally-safe-$k-steps.mat"
    matwrite(filename, Dict("indVmin" => indVmin, "indVmax" => indVmax, "policy" => policy))
end