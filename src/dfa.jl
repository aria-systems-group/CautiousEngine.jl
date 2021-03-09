struct DFA 
    states
    aps
    transitions
    accepting_state
    sink_state
    initial_state
end

"""
Transition function for any DFA.
"""
function δ(transitions, q, label)
    for relation in transitions 
        if relation[1] == q && (relation[2] == label || relation[2] == "true")
            return relation[3]
        end
    end
    return -1
end 

"""
Calculate Transition reward
"""
function distance_from_accept_state(dfa, dfa_state)

    @assert dfa_state ∈ dfa.states

    if dfa_state == dfa.sink_state
        return Inf
    end

    if dfa_state == dfa.accepting_state
        return 0
    end

    states_to_try = [dfa.accepting_state]
    dist = 1

    while !isempty(states_to_try)
        for state in states_to_try
            for relation in dfa.transitions
                if relation[3] == state 
                    if relation[1] == dfa_state 
                        # @info relation[1], state
                        return dist
                    else
                        push!(states_to_try, relation[1])
                        # @info states_to_try
                    end
                end
            end
            dist += 1
            setdiff!(states_to_try, state)
        end
    end

    return dist
end