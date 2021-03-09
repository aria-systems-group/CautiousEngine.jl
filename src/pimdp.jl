mutable struct PIMDP
    states
    actions
    Pbounds
    imdp_label_dict
    accepting_labels
    sink_labels
    initial_state
    policy
    state_history
    trajectory 
end

" 
Construct a PIMDP from an MDP by directly specifying unsafe states.
"
# TODO: Need multiple ways to do this. This is a hackish way for now.
function construct_PIMDP_from_IMDP(imdp, unsafe_states)

    # Problem 1: What are the accepting states?
    # Answer: the accepting states are the unsafe states. Then, we perform verification over the unsafe states and take the complementary
    # TODO: This is very hackkkkky
    sink_labels = zeros(1, length(imdp.states))
    sink_labels[unsafe_states] .= 1
    pimdp = PIMDP(imdp.states, imdp.actions, imdp.Pbounds, nothing, sink_labels, [],
                  nothing, nothing, nothing, nothing)

    return pimdp
end

"""
Construct a product IMDP from a (multimode) IMDP and DFA.
"""
function construct_DFA_IMDP_product(dfa, imdp)
    # Intial state
    qinit = δ(dfa.transitions, dfa.initial_state, imdp.labels[imdp.initial_state])
    sizeQ = length(dfa.states)
    pinit = (imdp.initial_state, qinit)
    dfa_acc_state = dfa.accepting_state
    dfa_sink_state = dfa.sink_state

    # Transition matrix size will be Nx|Q| -- or the original transition matrix permutated with the number of states in DFA
    Pbounds = imdp.Pbounds
    N, M = size(Pbounds)
    Pbounds_new = hcat([[0.0..0.0 for i=1:N*sizeQ] for j=1:M*sizeQ]...)

    pimdp_states = []
    pimdp_actions = imdp.actions
    sizeA = length(pimdp_actions)
    @debug "Actions: " pimdp_actionsA
    for s in imdp.states
        for q in dfa.states
            new_state = (s, q)
            push!(pimdp_states, new_state)
            # Check for transitions to other states
        end
    end

    for sq in pimdp_states
        for a in pimdp_actions
            for sqp in pimdp_states
                qp_test = δ(dfa.transitions, sq[2], imdp.labels[sq[1]])
                if qp_test == sqp[2]
                    # Get the corresponding entry of the Pbounds matrix
                    row_idx = (sq[1]-1)*sizeQ*sizeA + (sq[2]-1)*sizeA + a 
                    if (sq[2] == dfa_acc_state && sqp[2] == dfa_acc_state) || (sq[2] == dfa_sink_state && sqp[2] == dfa_sink_state)
                        # Flush out the old probabilities
                        col_idx = (sq[1]-1)*sizeQ + sq[2]  
                        Pbounds_new[row_idx, :] .= [0.0..0.0]
                        Pbounds_new[row_idx, col_idx] = 1.0..1.0
                    else
                        col_idx = (sqp[1]-1)*sizeQ + sqp[2]
                        Pbounds_new[row_idx, col_idx] = Pbounds[(sq[1]-1)*sizeA + a, sqp[1]]
                    end
                end
            end
        end
    end
   
    labels = zeros(1, M*sizeQ)
    labels[dfa_acc_state:sizeQ:M*sizeQ] .= 1

    sink_labels = zeros(1, M*sizeQ)
    sink_labels[dfa_sink_state:sizeQ:M*sizeQ] .= 1

    pimdp = PIMDP(pimdp_states, imdp.actions, Pbounds_new, imdp.labels, labels, sink_labels, nothing, nothing, nothing, nothing)
    return pimdp 
end

"""
Write a PIMDP object to file in a format consistent with Morteza's synthesis tool.
"""
function write_pimdp_to_file(pimdp, filename)
    open(filename, "w") do f
        state_num = length(pimdp.states)
        action_num =length(pimdp.actions)
        @printf(f, "%d \n", state_num)
        @debug "Length actions: " length(pimdp.actions)
        @printf(f, "%d \n", length(pimdp.actions))
        # Get number of accepting states from the labels vector
        @printf(f, "%d \n", sum(pimdp.accepting_labels))
        pimdp_acc_labels = pimdp.accepting_labels
        acc_states = findall(>(0.), pimdp_acc_labels[:])
        [@printf(f, "%d ", acc_state-1) for acc_state in acc_states]
        @printf(f, "\n")
        pimdp_sink_labels = pimdp.sink_labels
        sink_states = findall(>(0.), pimdp_sink_labels[:])

        for i=1:state_num
            if isnothing(sink_states) || !(i∈sink_states)
                for action in pimdp.actions
                    row_idx = (i-1)*action_num + action
                    ij = findall(>(0.), maximum.(pimdp.Pbounds[row_idx, :]))   
                    # Something about if the upper bound is less than one? Perhaps for numerical issues?
                    @debug action, i
                    psum = sum(maximum.(pimdp.Pbounds[row_idx, :]))
                    # @assert sum(maximum.(pimdp.Pbounds[row_idx, :])) >= 1
                    psum >= 1 ? nothing : throw(AssertionError("Bad max sum: $psum")) 
                    for j=ij
                        @printf(f, "%d %d %d %f %f", i-1, action-1, j-1, minimum(pimdp.Pbounds[row_idx, j]), maximum(pimdp.Pbounds[row_idx, j]))
                        if (i < state_num || j < ij[end] || action < action_num)
                            @printf(f, "\n")
                        end
                    end
                end
            else
                [@printf(f, "%d %d %d %f %f\n", i-1, j, i-1, 1.0, 1.0) for j=0:action_num-1]
            end
        end
    end
end

"""
Read a PIMDP from file.
"""
function read_pimdp_from_file(filename)
    pimdp = nothing
    open(filename, "r") do f
        num_states = tryparse(Int, readline(f))
        num_actions = tryparse(Int, readline(f)) 
        num_sink_states = tryparse(Int, readline(f))
        accept_states = parse.(Int, split(readline(f))) .+ 1
        Pbounds = hcat([[0.0..0.0 for i=1:num_states*num_actions] for j=1:num_states]...) 
        while !eof(f)
            row_split = split(readline(f))
            if length(row_split) > 0
                q0 = parse(Int, row_split[1])
                a = parse(Int, row_split[2]) + 1
                qt = parse(Int, row_split[3])
                plow = parse(Float64, row_split[4])
                phigh = parse(Float64, row_split[5])
                row_idx = q0*num_actions + a 
                col_idx = qt + 1
                Pbounds[row_idx, col_idx] = plow..phigh
            end
        end

        pimdp = PIMDP(collect(1:num_states), collect(1:num_actions), Pbounds, 
        nothing, accept_states, [], nothing, nothing, nothing, nothing)
    end

    return pimdp
end

"""
Quickly Validate a (P)IMDP object.
"""
function validate_pimdp(pimdp)
    for row in eachrow(pimdp.Pbounds)
        @assert sum(maximum.(row)) >= 1
        @assert sum(minimum.(row)) <= 1
    end
end


"""
Determine the IMDP state from the current PIMDP index and # of DFA states
"""
function pimdp_col_to_imdp_state(col, num_dfa_states)
    return ceil(col/num_dfa_states)
end


"""
% computes the true distribution given a range of distributions
    
% input:    pmin: lower bound transition probabilities
%           pmax: upper bound transition probabilities
%           indSorted: index of states according to their sorted list
% output:
%           p: true transition probability
"""
function get_true_transition_probabilities(pmin, pmax, indSorted)
    p = zeros(1,length(pmax))
    used = sum(pmin)
    remain = 1 - used
    
    for i in indSorted
        if pmax[i] <= (remain + pmin[i])
            p[i] = pmax[i]
        else
            p[i] = pmin[i] + remain
        end
        remain = maximum([0., remain - (pmax[i] - pmin[i]) ])   
    end

    return p
end