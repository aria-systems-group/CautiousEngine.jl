struct IMDP
    states
    actions
    Pbounds
    labels
    sink_labels
    initial_state
end

using IntervalSets

"""
Create an IMDP from a list of result directories.
"""
function create_imdp_from_result_dirs(res_dirs, exp_dir)
    sys_dir, minPr, maxPr = create_switched_system(res_dirs, exp_dir)
    N, M = size(minPr)
    Pbounds = hcat([[minPr[i, j]..maxPr[i,j] for i=1:N] for j=1:M]...)
    Q = collect(1:M)
    A = collect(1:Int(N/M))
    @debug "Actions: " A
    labels = Dict()
    imdp = IMDP(Q, A, Pbounds, labels, nothing, 1)
    return imdp
end

"
Create an IMDP from a set of transition bound matrices.
"
function create_single_mode_imdp(minPr, maxPr)
    N, _ = size(minPr)
    Pbounds = hcat([[minPr[i, j]..maxPr[i,j] for i=1:N] for j=1:N]...)
    Q = collect(1:N)
    A = [1]
    labels = Dict()
    imdp = IMDP(Q, A, Pbounds, labels, nothing, 1)
    return imdp 
end

function create_imdp_labels(labels_fn, imdp, extent_file)
    res = BSON.load(extent_file)
    imdp_state_extents = res[:extents]
    for state_key in keys(imdp_state_extents)
        state = imdp_state_extents[state_key]
        if state_key == -11
            imdp.labels[length(imdp_state_extents)] = labels_fn(state, unsafe=true) 
        else
            imdp.labels[state_key] = labels_fn(state)
        end
    end
end