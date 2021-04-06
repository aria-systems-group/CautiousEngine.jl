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
function create_simple_imdp(minPr, maxPr)
    N, _ = size(minPr)
    Pbounds = hcat([[minPr[i, j]..maxPr[i,j] for i=1:N] for j=1:N]...)
    Q = collect(1:N)
    A = [1]
    labels = Dict()
    for i=1:N-1
        labels[i] = "safe"
    end
    labels[N] = "!safe"
    imdp = IMDP(Q, A, Pbounds, labels, nothing, 1)
    return imdp 
end

function create_imdp_labels(labels_fn, imdp, extent_file::String)
    res = BSON.load(extent_file)
    imdp_state_extents = res[:extents]
    create_imdp_labels(labels_fn, imdp, imdp_state_extents)
end

function create_imdp_labels(labels_fn, imdp, imdp_state_extents::Dict)
    for state_key in keys(imdp_state_extents)
        state = imdp_state_extents[state_key]
        if state_key == -11 || state_key == length(keys(imdp_state_extents))
            imdp.labels[length(imdp_state_extents)] = labels_fn(state, unsafe=true) 
        else
            imdp.labels[state_key] = labels_fn(state)
        end
    end
end

function write_imdp_to_file_bounded(imdp, Qyes, Qno, filename)
    open(filename, "w") do f
        state_num = length(imdp.states)
        action_num =length(imdp.actions)
        @printf(f, "%d \n", state_num)
        @debug "Length actions: " length(imdp.actions)
        @printf(f, "%d \n", length(imdp.actions))
        # Get number of accepting states from the labels vector
        acc_states = Qyes 
        @printf(f, "%d \n", length(acc_states))
        [@printf(f, "%d ", acc_state-1) for acc_state in acc_states]
        @printf(f, "\n")
        sink_states = Qno 

        for i=1:state_num
            if isnothing(sink_states) || !(i∈sink_states)
                for action in imdp.actions
                    row_idx = (i-1)*action_num + action
                    ij = findall(>(0.), maximum.(imdp.Pbounds[(i-1)*action_num + action, :]))   
                    # Something about if the upper bound is less than one? Perhaps for numerical issues?
                    @debug action, i
                    psum = sum(maximum.(imdp.Pbounds[row_idx, :]))
                    psum >= 1 ? nothing : throw(AssertionError("Bad max sum: $psum")) 
                    for j=ij
                        @printf(f, "%d %d %d %f %f", i-1, action-1, j-1, minimum(imdp.Pbounds[row_idx, j]), maximum(imdp.Pbounds[row_idx, j]))
                        if (i < state_num || j < ij[end] || action < action_num)
                            @printf(f, "\n")
                        end
                    end
                end
            else
                @printf(f, "%d %d %d %f %f", i-1, 0, i-1, 1.0, 1.0)
                if i<state_num
                    @printf(f, "\n")
                end
            end
        end
    end
end

function create_dot_graph(imdp::IMDP, filename::String)
    open(filename, "w") do f
        println(f, "digraph G {")
        println(f, "  rankdir=LR\n  node [shape=\"circle\"]\n  fontname=\"Lato\"\n  node [fontname=\"Lato\"]\n  edge [fontname=\"Lato\"]")
        println(f, "  size=\"8.2,8.2\" node[style=filled,fillcolor=\"#FDEDD3\"] edge[arrowhead=vee, arrowsize=.7]")

        # Initial State
        @printf(f, "  I [label=\"\", style=invis, width=0]\n  I -> %d\n", imdp.initial_state)

        for state in imdp.states
            @printf(f, "  %d [label=<q<SUB>%d</SUB>>, xlabel=<%s>]\n", state, state, imdp.labels[state])
            
            for action in imdp.actions
                row_idx = (state-1)*length(imdp.actions) + action

                for idx in findall(>(0.), maximum.(imdp.Pbounds[row_idx, :]))
                    state_p = idx
                    @printf(f, "  %d -> %d [label=<a<SUB>%d</SUB>: %.1f-%.1f >]\n", state, state_p, action, minimum(imdp.Pbounds[row_idx,state_p]), maximum(imdp.Pbounds[row_idx,state_p]))
                end
            end
        end
        println(f, "}")
    end
end

"""
Get IMDP Ring states for a gridworld
"""
function get_imdp_ring_states(imdp_state_idx, num_stats_per_col, num_imdp_states; col_ordering=true)
    # TODO: enfore column ordering 
    ring_imdp_states = [imdp_state_idx + 1, imdp_state_idx - 1,
                        imdp_state_idx - num_stats_per_col, imdp_state_idx - num_stats_per_col + 1, imdp_state_idx - num_stats_per_col - 1,
                        imdp_state_idx + num_stats_per_col, imdp_state_idx + num_stats_per_col + 1, imdp_state_idx + num_stats_per_col - 1]

   # To do: Better filter

   # Remove the negative states 
   ring_imdp_states = ring_imdp_states[ring_imdp_states .> 1]
   ring_imdp_states = ring_imdp_states[ring_imdp_states .< num_imdp_states]
   return ring_imdp_states
end