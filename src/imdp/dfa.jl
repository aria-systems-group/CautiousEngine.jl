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

    # TODO: If no transitions, stay in the current state - this is a workaround to use Spot.jl
    return q
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

""" 
Construct DFA object using Spot
"""
function construct_DFA_from_LTL(ltl_form)
    # ltl_form = ltl"$ltl_string"
    translator = LTLTranslator(deterministic=true)
    spot_automaton = Spot.translate(translator, ltl_form)

    transitions = []
    for (edge, label) in zip(get_edges(spot_automaton), get_labels(spot_automaton))
        label = string(label) == "1" ? "true" :  create_full_label(atomic_propositions(spot_automaton), label)
        new_transition = (edge[1], label, edge[2])
        push!(transitions, new_transition)
    end

    # TODO: add transition to the/a violating state
    # ! It doesn't exist here currently
    accept_states = get_rabin_acceptance(spot_automaton)
    accept_state = collect(accept_states[1][2])[1]

    dfa = DFA(collect(1:num_states(spot_automaton)), string.(atomic_propositions(spot_automaton)), transitions, accept_state, nothing, get_init_state_number(spot_automaton))

    return dfa
end

function create_full_label(props, label)
    # all props
    label_props = atomic_prop_collect(label)
    if props == label_props
        return string(label)
    end

    # rebuild
    label = string(label)
    new_label = ""
    for prop in string.(props)
        if occursin("!$prop", label)
            new_label = "$new_label & !$prop"
        elseif occursin("$prop", label)
            new_label = "$new_label & $prop"
        else
            new_label = "$new_label & !$prop"
        end
    end
    new_label = new_label[4:end]
    return new_label
end
