using Printf 

function create_graph_from_imdp(imdp, filename)
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

function create_graph_from_dfa(dfa, filename)
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

function create_graph_from_pimdp(pimdp, filename)
    open(filename, "w") do f
        println(f, "digraph G {")
        println(f, "  rankdir=LR\n  node [shape=\"circle\"]\n  fontname=\"Lato\"\n  node [fontname=\"Lato\"]\n  edge [fontname=\"Lato\"]")
        println(f, "  size=\"8.2,8.2\" node[style=filled,fillcolor=\"#FDEDD3\"] edge[arrowhead=vee, arrowsize=.7]")

        # Initial State
        # TODO: Don't hardcode intial state
        @printf(f, "  I [label=\"\", style=invis, width=0]\n  I -> %d\n", 1)

        i = 1
        for state in pimdp.states
            @printf(f, "  %d [label=<(s<SUB>%d</SUB>,q<SUB>%d</SUB>)>, xlabel=<%d>]\n", i, state[1], state[2], pimdp.accepting_labels[(state[1]-1)*3 + state[2]])
            
            for action in pimdp.actions
                q_size = 3
                row_idx = (state[1]-1)*length(pimdp.actions)*q_size + (state[2]-1)*length(pimdp.actions) + action

                for idx in findall(>(0.), maximum.(pimdp.Pbounds[row_idx, :]))
                    state_p = pimdp.states[idx]
                    col_idx = (state_p[1]-1)*q_size + state_p[2]
                    # TODO: fix the second index to be correct
                    @printf(f, "  %d -> %d [label=<a<SUB>%d</SUB>: %.1f-%.1f >]\n", i, idx, action, minimum(pimdp.Pbounds[row_idx,col_idx]), maximum(pimdp.Pbounds[row_idx,col_idx]))
                end
            end
            i+=1
        end
        println(f, "}")
    end
end
