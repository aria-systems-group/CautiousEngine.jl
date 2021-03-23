using CautiousEngine

@testset "dfa function tests" begin
    # Simple reach-avoid Specification
    dfa_states = collect(1:3)
    dfa_props = ["a", "b"]
    dfa_transitions = [(1, "!a∧!b", 1),
                        (1, "a∧!b", 2),
                        (2, "true", 2),
                        (1, "!a∧b", 3),
                        (3, "true", 3)]
    dfa_accepting_state = 2
    dfa_sink_state = 3
    dfa_initial_state = 1

    dfa = CautiousEngine.DFA(dfa_states, dfa_props, dfa_transitions, dfa_accepting_state, dfa_sink_state, dfa_initial_state)

    @test CautiousEngine.get_next_dfa_labels(dfa, 1) == ["a∧!b"]
    @test CautiousEngine.get_next_dfa_labels(dfa, 2) == []
    @test CautiousEngine.get_next_dfa_labels(dfa, 3) == []

end