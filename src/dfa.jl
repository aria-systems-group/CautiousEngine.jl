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
