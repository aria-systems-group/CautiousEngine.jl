function simulate_system(x0, modes, steps, imdp, pimdp, dfa, extents, policy; noise_dist=nothing)
    reset_pimdp(x0, imdp, dfa, pimdp, extents, policy)
    for i=1:steps
        check_pimdp_exit_conditions(pimdp) ? break : nothing
        x_c = pimdp.trajectory[end]
        mode = pimdp.policy[map_pimdp_state_to_index(pimdp, pimdp.state_history[end])]
        x_n = modes[Int(mode)](x_c)
        if !isnothing(noise_dist)
            x_n += rand(noise_dist, size(x_n))
        end
        propogate_pimdp_trajectory(pimdp, dfa, extents, x_n)
    end
    return 
end

function check_pimdp_exit_conditions(pimdp)
    current_state = pimdp.state_history[end]
    current_state_idx = map_pimdp_state_to_index(pimdp, current_state)

    if current_state_idx ∈ findall(>(0.), pimdp.accepting_labels[:])
        @info "We have reached a positive accepting state. Exiting."
        return true
    elseif current_state_idx ∈ findall(>(0.), pimdp.sink_labels[:])
        @info "We have reached a negative sink state. Booo. Exiting."
        return true
    end
    return false 
end

function map_pimdp_state_to_index(pimdp, state)
    index = nothing
    for (i, test_state) in enumerate(pimdp.states)
        test_state == state ? index = i : nothing 
        !isnothing(index) ? break : nothing
    end
    @assert !isnothing(index)
    return index
end

function propogate_pimdp_trajectory(pimdp, dfa, extent_dict, x_new)
    new_extent = nothing
    for extent_id in keys(extent_dict)
        extent_id == -11 ? continue : nothing
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

function reset_pimdp(x0, imdp, dfa, pimdp, extent_dict, policy)
    init_extent = nothing
    for extent_id in keys(extent_dict)
        extent_id == -11 ? continue : nothing
        extent = extent_dict[extent_id]
        if sum([extent[dim][1]<=x0[i]<=extent[dim][2] for (i, dim) in enumerate(keys(extent))]) == length(x0)
            init_extent = extent_id 
            @info init_extent
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