function simulate_system(x0, modes, steps, imdp, pimdp, dfa, extents, policy)
    x0 = reshape(x0, length(x0), 1)
    reset_pimdp(x0, imdp, dfa, pimdp, extents, policy)
    for i=1:steps
        check_pimdp_exit_conditions(pimdp) ? break : nothing
        x_c = pimdp.trajectory[end]
        mode = pimdp.policy[map_pimdp_state_to_index(pimdp, pimdp.state_history[end])]
        x_n = modes[Int(mode)](x_c)
        propogate_pimdp_trajectory(pimdp, dfa, extents, x_n)
    end
    return 
end

"""
Check the PIMDP exit conditions.
"""
function check_pimdp_exit_conditions(pimdp::PIMDP)
    current_state = pimdp.state_history[end]
    current_state_idx = map_pimdp_state_to_index(pimdp, current_state)

    if current_state_idx ∈ findall(>(0.), pimdp.accepting_labels[:])
        @debug "We have reached a positive accepting state. Exiting."
        return true
    elseif current_state_idx ∈ findall(>(0.), pimdp.sink_labels[:])
        @debug "We have reached a negative sink state. Booo. Exiting."
        return true
    end
    return false 
end

"""
Check the PIMDP exit conditions.
"""
function check_pimdp_exit_conditions(pimdp::PIMDP, num_dfa_states::Int)
    current_state = pimdp.state_history[end]
    current_state_idx = map_pimdp_state_to_index(current_state, num_dfa_states)

    acc_states = findall(>(0.), pimdp.accepting_labels[:])
    sink_states = findall(>(0.), pimdp.sink_labels[:]) 
    if current_state_idx ∈ acc_states 
        @debug "We have reached a positive accepting state. Exiting."
        return 1 
    elseif current_state_idx ∈ sink_states
        @debug "We have reached a negative sink state. Booo. Exiting."
        return 0 
    end
    return nothing 
end

"""
Find the enumeration of the state tuple in terms of the PIMDP state.
"""
function map_pimdp_state_to_index(pimdp::PIMDP, state::Tuple)
    index = nothing
    for (i, test_state) in enumerate(pimdp.states)
        test_state == state ? index = i : nothing 
        !isnothing(index) ? break : nothing
    end
    @assert !isnothing(index)
    return index
end

"""
Find the enumeration of the state tuple in terms of the PIMDP state. 
"""
function map_pimdp_state_to_index(state::Tuple, num_dfa_states::Int)
   index = (state[1]*num_dfa_states) - num_dfa_states + state[2]
   return index 
end

"""
Find the row of the Pbounds matrix given the PIMDP state and action #
"""
function map_pimdp_state_to_Pmat_row(state::Tuple, num_dfa_states::Int, num_actions::Int, action::Int)
    row_idx = (state[1]-1)*num_dfa_states*num_actions + (state[2]-1)*num_actions + action 
    return row_idx
end

function propogate_pimdp_trajectory(pimdp, dfa, extent_dict, x_new)
    new_extent = nothing
    for extent_id in keys(extent_dict)
        extent_id == extent_dict.count ? continue : nothing
        extent = extent_dict[extent_id]
        if sum([extent[dim][1]<=x_new[i]<=extent[dim][2] for (i, dim) in enumerate(keys(extent))]) == length(x_new)
            new_extent = extent_id 
            break
        end
    end

    if isnothing(new_extent)
        qnew = dfa.sink_state
        new_extent = length(keys(extent_dict)) 
    else
        qnew = δ(dfa.transitions, pimdp.state_history[end][2], pimdp.imdp_label_dict[new_extent])
    end
    push!(pimdp.state_history, (new_extent, qnew))
    push!(pimdp.trajectory, x_new)
end

"""
Get IMDP State ID from system state
"""
function get_imdp_state_id(extent_dict, x_new)
    new_extent = nothing
    for extent_id in keys(extent_dict)
        extent_id == extent_dict.count ? continue : nothing
        extent = extent_dict[extent_id]
        if sum([extent[dim][1]<=x_new[i]<=extent[dim][2] for (i, dim) in enumerate(keys(extent))]) == length(x_new)
            new_extent = extent_id 
            break
        end
    end

    if isnothing(new_extent)
        new_extent = length(keys(extent_dict))
    end
    return new_extent
end

"""
Test PIMDP Propogation without saving
"""
function propogate_pimdp_trajectory_test(pimdp, dfa, extent_dict, x_new)
    new_extent = nothing
    for extent_id in keys(extent_dict)
        extent_id == extent_dict.count ? continue : nothing
        extent = extent_dict[extent_id]
        if sum([extent[dim][1]<=x_new[i]<=extent[dim][2] for (i, dim) in enumerate(keys(extent))]) == length(x_new)
            new_extent = extent_id 
            break
        end
    end

    if isnothing(new_extent)
        qnew = dfa.sink_state
        new_extent = length(keys(extent_dict)) 
    else
        qnew = δ(dfa.transitions, pimdp.state_history[end][2], pimdp.imdp_label_dict[new_extent])
    end

    # reward = 0.
    # if qnew == dfa.sink_state
    #     reward = -Inf
    # elseif qnew == dfa.accepting_state
    #     reward = Inf
    # elseif qnew > pimdp.state_history[end][2]
    #     reward = 1000
    # end
    
    # return reward
    return (new_extent, qnew)
end

"""
Reward based on the worst-possible next DFA state
"""
function calculate_worst_dfa_reward(pimdp, dfa, possible_extents)
    best_reward = Inf
    for extent in possible_extents
        qnew = δ(dfa.transitions, pimdp.state_history[end][2], pimdp.imdp_label_dict[extent])
        reward = -distance_from_accept_state(dfa, qnew)
        if reward < best_reward
            best_reward = reward
        end
    end

    return best_reward
end

"""
Reward based on the best-possible next DFA state
"""
function calculate_best_dfa_reward(pimdp, dfa, possible_extents)
    best_reward = -Inf
    for extent in possible_extents
        qnew = δ(dfa.transitions, pimdp.state_history[end][2], pimdp.imdp_label_dict[extent])
        reward = -distance_from_accept_state(dfa, qnew)
        if reward > best_reward
            best_reward = reward
        end
    end

    return best_reward
end

"""
Reward based on smallest distance to accept state in PIMDP
"""
function reset_pimdp(x0, imdp, dfa, pimdp, extent_dict, policy)
    init_extent = nothing
    for extent_id in keys(extent_dict)
        extent_id == extent_dict.count ? continue : nothing
        extent = extent_dict[extent_id]
        if sum([extent[dim][1]<=x0[i]<=extent[dim][2] for (i, dim) in enumerate(keys(extent))]) == length(x0)
            init_extent = extent_id 
            break
        end
    end
    @assert !isnothing(init_extent)

    qinit = δ(dfa.transitions, dfa.initial_state, imdp.labels[init_extent])
    initial_pimdp_state = (init_extent, qinit)
    pimdp.initial_state = initial_pimdp_state
    pimdp.state_history = [initial_pimdp_state]
    pimdp.trajectory = [x0]
    pimdp.policy = policy
end