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

function create_dot_graph(dfa::DFA, filename::String)
    open(filename, "w") do f
        println(f, "digraph G {")
        println(f, "  rankdir=LR\n  node [shape=\"circle\"]\n  fontname=\"Lato\"\n  node [fontname=\"Lato\"]\n  edge [fontname=\"Lato\"]")
        println(f, "  size=\"8.2,8.2\" node[style=filled,fillcolor=\"#FDEDD3\"] edge[arrowhead=vee, arrowsize=.7]")

        # Initial State
        @printf(f, "  I [label=\"\", style=invis, width=0]\n  I -> %d\n", dfa.initial_state)

        for state in dfa.states
            if state == dfa.accepting_state
                @printf(f, "  %d [shape=\"doublecircle\", label=<z<SUB>%d</SUB>>]\n", state, state)
                @printf(f, "  %d -> %d [label=<true>]\n", state, state)
            else
                @printf(f, "  %d [label=<z<SUB>%d</SUB>>]\n", state, state)
                for transition in dfa.transitions
                    if transition[1] == state
                        # Assume 1 label for now
                        @printf(f, "  %d -> %d [label=<%s>]\n", state, transition[3], transition[2])
                    end
                end

            end

        end
        println(f, "}")
    end
end

"""
Given the current dfa state, return labels that will induce a positive transition
"""
function get_next_dfa_labels(dfa, dfa_state)
    labels = []
    for relation in dfa.transitions
        if relation[1] == dfa_state &&  relation[1] != relation[3] && relation[3] != dfa.sink_state
            push!(labels, relation[2])
        end
    end
    return labels
end

"""
Get all labels that transition to sink state
"""
function get_dfa_sink_labels(dfa)
    labels = []
    for relation in dfa.transitions
        if relation[1] != relation[3] && relation[3] == dfa.sink_state
            push!(labels, relation[2])
        end
    end
    return labels
end

"""
Check if current state can transition to accept state
"""
function possible_accept_transition(dfa, dfa_state, symbol)
    for relation in dfa.transitions
        if dfa_state == relation[1] && symbol == relation[2] && relation[3] == dfa.accepting_state
            return true 
        end
    end
    return false
end