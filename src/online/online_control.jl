# Methods related to online control of dynamic systems

"""
Objects relevant to online control problem
"""
struct OnlineControlProblem
    initial_state
    control_modes
    # TODO: We don't need this IMDP here, get rid of it.
    imdp::IMDP
    pimdp::PIMDP
    dfa::DFA
    extents
    post_image_data
    control_policy
    local_gps::Dict
    global_gps::Dict # stores the paths to the full GPs generated offline
    reward_type::Tuple
    # TODO: There is probably a better place for this somewhere else! 
    primary_reward_history
    secondary_reward_history
    online_data # Stores the data collected online
    global_gp_infos # Stores the persistent info dicts
end

"""
Setup the online control problem.
"""
function setup_online_control_problem(x0, modes, imdp, pimdp, dfa, extents, posts, policy, local_gps, global_gps, reward_type)
    reset_pimdp(x0, imdp, dfa, pimdp, extents, policy)

    # Setup local gp data
    online_data = Dict()
    for i=1:length(modes)
        online_data[i] = Dict("x_online" => [], "y_online" => Dict(i=>[] for i in keys(extents[1])))
    end
    # online_data["x_online"] = []
    # online_data["y_online"] = Dict(i=>[] for i in keys(extents[1]))

    global_gp_infos = Dict()

    return OnlineControlProblem(x0, modes, imdp, pimdp, dfa, extents, posts, policy, local_gps, global_gps, reward_type, [], [], online_data, global_gp_infos)
end

"""
One step of a control loop with an offline policy.
"""
function online_control_loop(ocp::OnlineControlProblem)
    res = check_pimdp_exit_conditions(ocp.pimdp, length(ocp.dfa.states)) 
    if !isnothing(res) 
        return res
    end
    x_c = ocp.pimdp.trajectory[end]
    mode = ocp.pimdp.policy[map_pimdp_state_to_index(ocp.pimdp.state_history[end], length(ocp.dfa.states) )]
    x_n = ocp.control_modes[Int(mode)](x_c)
    propogate_pimdp_trajectory(ocp.pimdp, ocp.dfa, ocp.extents, x_n)
    return nothing
end

"""
Concatonate Online and Offline Data
"""
function concatonate_data_points(x_data_offline, y_data_offline, x_data_online, y_data_online)
    x_cat = [x_data_offline x_data_online...]
    y_cat = [y_data_offline; y_data_online]
    return x_cat, y_cat
end

using NearestNeighbors
"""
Create local GP for an extent.
"""
function get_local_data_knn(center, x_data, y_data; num_neighbors = 50)
    num_neighbors = minimum([num_neighbors, size(x_data,2)])
    kdtree = KDTree(x_data);
    sub_idx, _ = knn(kdtree, center, num_neighbors, true)

    # WHY THE HELL IS IT DOING THIS
    if typeof(sub_idx) == Array{Array{Int64,1},1}
        sub_idx = sub_idx[1]
    end

    sub_x_data = x_data[:, sub_idx]
    sub_y_data = y_data[sub_idx]
    return sub_x_data, sub_y_data 
end

"""
One step optimal control loop with an offline policy and GP update.
"""
function online_control_loop_optimal(ocp::OnlineControlProblem; 
                                     update_gp_flag=false, gp_use_key="global", update_pimdp_ring_states=false,
                                     pimdp_reward_dict=nothing)
    etime = @elapsed begin
        # Check if an exit condition has been reached
        res = check_pimdp_exit_conditions(ocp.pimdp, length(ocp.dfa.states)) 
        if !isnothing(res) 
                
            return res
        end

        x_c = ocp.pimdp.trajectory[end]     # Current real state
        q_c = ocp.pimdp.state_history[end]
        pimdp_state = map_pimdp_state_to_index(q_c, length(ocp.dfa.states))  # Current PIMDP state index      

        best_mode, best_reward, gp_info_dict = get_best_action(ocp, gp_use_key=gp_use_key, pimdp_reward_dict=pimdp_reward_dict)
        x_n = ocp.control_modes[best_mode](x_c)

        propogate_pimdp_trajectory(ocp.pimdp, ocp.dfa, ocp.extents, x_n)

        @debug "BEST ACTION: " best_mode
        @debug "BEST REWARD: " best_reward

        push!(ocp.primary_reward_history, best_reward[1])
        push!(ocp.secondary_reward_history, best_reward[2])

        # TODO: Move this to its own fcn / subfcns
        if update_gp_flag && isnothing(check_pimdp_exit_conditions(ocp.pimdp, length(ocp.dfa.states)))
            dim_keys = keys(gp_info_dict)
            push!(ocp.online_data[best_mode]["x_online"], x_c)
            # Update the GP with a new datapoint
            for (i, dim_key) in enumerate(dim_keys)
                y_meas = (x_n[i] - x_c[i])
                push!(ocp.online_data[best_mode]["y_online"][dim_key], y_meas)
                add_datapoint_to_gp(gp_info_dict[dim_key].gp, (x_c, y_meas))
                update_gp_info(gp_info_dict[dim_key])
            end

            pimdp_states_to_update = [q_c, ]
            if update_pimdp_ring_states
                new_states = get_pimdp_ring_states(q_c)

                for new_state in new_states
                    push!(pimdp_states_to_update, new_state)
                end
            end

            for pimdp_state in pimdp_states_to_update

                # Update transition probability to some states
                extent_id = pimdp_state[1]
                regions_to_update = []
                # pimdp_row = map_pimdp_state_to_Pmat_row(q_c, length(ocp.dfa.states), , best_mode)
                pimdp_rows = get_pimdp_row_idx_same_imdp(pimdp_state, ocp.dfa, length(ocp.pimdp.actions), best_mode)
                @assert length(pimdp_rows) > 0 

                # og_post = ocp.post_image_data[best_mode][pimdp_state[1]]

                # Get the new post-image of the GP
                current_extent = ocp.extents[extent_id]
                lb = [current_extent[dim_key][1] for dim_key in dim_keys]
                ub = [current_extent[dim_key][2] for dim_key in dim_keys]

                ## TODO: THIS IS MESSY
                region_post_dict_new = bound_extent(current_extent, lb, ub, gp_info_dict, dim_keys, Dict("x1" => [1., 1.], "x2" => [1., 1.]); known_part_flag=true)
                dummy_region_post_dict = Dict(extent_id => region_post_dict_new)
                
                for pimdp_row_idx in pimdp_rows
                    prow_min = ocp.pimdp.Pmin[pimdp_row_idx, :] 
                    prow_max = ocp.pimdp.Pmax[pimdp_row_idx, :]

                    # Find all states with interval bounds
                    regions_to_update = findall(x->x>0, prow_min .!= prow_max)
                    pimdp_row_vals_min = prow_min[regions_to_update]
                    pimdp_row_vals_max = prow_max[regions_to_update]

                    noise_proc = Truncated(Normal(0, 0.01), -0.01, 0.01)
                    η = 0.01
                    ϵ = -1

                    for target in regions_to_update
                        target_id = pimdp_col_to_imdp_state(target, length(ocp.dfa.states))
                        region_pair = (extent_id, target_id)
                        frame_row = region_pair_transitions(region_pair, ocp.extents, dummy_region_post_dict, ϵ, gp_info_dict, noise_proc, η)

                        if target_id == 1025.
                            minPrn = 1 - frame_row[4] 
                            maxPrn = 1 - frame_row[3]

                            if minPrn > prow_min[target]
                                prow_min[target] = minPrn
                            end
                            if maxPrn < prow_max[target]
                                prow_max[target] = maxPrn
                            end
                        else
                            if frame_row[4] < prow_max[target]
                                prow_max[target] = frame_row[4]
                            end
                            if frame_row[3] > prow_min[target]
                                prow_min[target] = frame_row[3]
                            end
                        end
                    end

                    ocp.pimdp.Pmin[pimdp_row_idx, :] = prow_min
                    ocp.pimdp.Pmax[pimdp_row_idx, :] = prow_max
                end
            end
        end
    end
    @debug "Total time: $etime seconds"
end

"""
Update PIMDP row
"""
function update_pimdp_transition_bound_row()
    nothing
end


"""
Get the relevant GP info
"""
function get_gp_info(ocp, mode, current_state; gp_use_case="global")
    # Begin all of the fancy local GP processing
    global_gp_key = "global-$mode"
            
    # gp_use_case = "use-local-fresh"
    if gp_use_case == "global"
        if global_gp_key in keys(ocp.global_gp_infos)
            global_gp_info_dict = ocp.global_gp_infos[global_gp_key]
        else
            gp_filename = ocp.global_gps[mode]
            gp_info = load_gp("$gp_filename")
            global_gp_info_dict = gp_info[:gp_info]  
            ocp.global_gp_infos[global_gp_key] = global_gp_info_dict
        end
        gp_info_dict = global_gp_info_dict 
    elseif gp_use_case == "local-persistent"
        nothing
        # gp_key = "$imdp_state-$best_mode"
        # if gp_key in keys(ocp.local_gps) 
        #     gp_info_dict = ocp.local_gps[gp_key]
        #     # gp_info_dict
        # else
        #     #gp_filename = filter(x -> occursin("-$imdp_state.bin", x), readdir(gp_dir))[1]
        #     # Load the global GP this could be optimized
        #     gp_filename = ocp.global_gps[best_mode]
        #     gp_info = load_gp("$gp_filename")
        #     # gp_info = load_gp("$gp_dir/$gp_filename")
        #     # gp_local = create_local_gp_knn(x_c, gp_info; num_neighbors = 100) 
        #     gp_info_dict = gp_info[:gp_info]  
        # end

    elseif gp_use_case == "local-fresh"
        if global_gp_key in keys(ocp.global_gp_infos)
            global_gp_info_dict = ocp.global_gp_infos[global_gp_key]
        else
            gp_filename = ocp.global_gps[mode]
            gp_info = load_gp("$gp_filename")
            global_gp_info_dict = gp_info[:gp_info]  
            ocp.global_gp_infos[global_gp_key] = global_gp_info_dict
        end

        # Create a fresh local GP from the global GP
        gp_info_dict = copy(global_gp_info_dict)
        x_offline = global_gp_info_dict["x1"].gp.x 
        x_online = ocp.online_data[mode]["x_online"]

        for dim_key in keys(gp_info_dict) 
            y_offline = global_gp_info_dict[dim_key].gp.y 
            y_online = ocp.online_data[mode]["y_online"][dim_key]        
            x_cat, y_cat = concatonate_data_points(x_offline, y_offline, x_online, y_online)
            x_nn, y_nn = get_local_data_knn(current_state, x_cat, y_cat, num_neighbors=75)

            gp_info_dict[dim_key].gp = train_gp_1dim(x_nn, y_nn; se_params=[0., 0.65])
            update_gp_info(gp_info_dict[dim_key])
        end
    else
        @error "Unknown GP use case: " + gp_use_case
    end

    return gp_info_dict
end

"""
Get the policy from a verification result.
"""
function get_policy_from_offline_synthesis(filename)
    policy = Dict()
    open(filename, "r") do f
        while !eof(f)
            row_split = split(readline(f))
            if length(row_split) > 0
                q0 = parse(Int, row_split[1]) + 1
                u = parse(Int, row_split[2]) + 1
                policy[q0] = u
            end
        end
    end
    return policy
end

"""
Gets the best action at each state.
"""
function get_best_action(ocp::OnlineControlProblem; gp_use_key="global", pimdp_reward_dict=nothing)
    # Get the current state and PIMDP state
    x_c = ocp.pimdp.trajectory[end]     # Current real state
    q_c = ocp.pimdp.state_history[end]  # Current IMDP (extent ID) / DFA state
    pimdp_state = map_pimdp_state_to_index(q_c, length(ocp.dfa.states))  # Current PIMDP state index      

    # Initialize the rewards and best control mode 
    best_reward = [-Inf, -Inf]
    best_ub = -Inf
    best_mode = nothing
    best_gp_info_dict = nothing

    # Check if the current policy gives us a shortcut
    if ocp.control_policy[pimdp_state, 2] == 1.0
        best_mode = Int(ocp.control_policy[pimdp_state, 1])
        best_reward[1] = 2
        best_gp_info_dict = get_gp_info(ocp, best_mode, x_c, gp_use_case=gp_use_key)
    else
        for mode = 1:length(ocp.control_modes)

            gp_info_dict = get_gp_info(ocp, mode, x_c, gp_use_case=gp_use_key)
            # Extract the relevant row of the PIMDP transition matrix
            P_row_idx = map_pimdp_state_to_Pmat_row(q_c, length(ocp.dfa.states), length(ocp.pimdp.actions), mode)
            prow_max = ocp.pimdp.Pmax[P_row_idx, :]

            # Find all states with interval bounds with non-zero UB
            regions_with_nonzero_ub = findall(x->x>0., prow_max) 
            @debug regions_with_nonzero_ub

            ##### Test out the point region transition fcn
            ϵ = -1.
            μ, σ = predict_one_step_full(gp_info_dict, x_c)
            # TODO: get rid of this noise dist implicit
            noise_sigma = 0.01
            dist = Truncated(Normal(0., noise_sigma), -noise_sigma, noise_sigma)

            non_trivial_regions = []
            prob_intervals = Dict()
            imdp_extent_ids = []
            nt_pimdp_states = []

            for reg in regions_with_nonzero_ub
                region_idx = pimdp_col_to_imdp_state(reg, length(ocp.dfa.states))
                push!(nt_pimdp_states, (region_idx, q_c[2]))
                push!(imdp_extent_ids, region_idx)
                res = point_region_transitions(μ, σ, region_idx, ocp.extents, ϵ, gp_info_dict; process_noise_dist=dist)
                if res[2] > 0.
                    push!(non_trivial_regions, map_pimdp_state_to_index((Int(region_idx), q_c[2]), length(ocp.dfa.states) ))
                    prob_intervals[map_pimdp_state_to_index((Int(region_idx), q_c[2]), length(ocp.dfa.states) )] = [res[1], res[2]]
                end
            end

            # Extract the potential actions based on the potential regions we could transition to.
            potential_actions = zeros(length(non_trivial_regions), 6) 
            total_reward = 0.
            
            max_reward = 0.
            for (i,region) in enumerate(non_trivial_regions)
                potential_actions[i, 1] = region
                potential_actions[i, 2:4] = ocp.control_policy[region,:]
                potential_actions[i, 5] = prob_intervals[region][1]
                potential_actions[i, 6] = prob_intervals[region][2]
                # total_reward += ocp.control_policy[region,2]
                # total_reward += ocp.control_policy[region,3]                    
            end
            # total_reward *= 1. / length(non_trivial_regions)

            # Calculate the worst expected lower bound of satisfaction
            sorted_actions = potential_actions[sortperm(potential_actions[:,3]), :] 
            region_idx = sorted_actions[:,1]
            minPrs = sorted_actions[:,5]
            maxPrs = sorted_actions[:,6]
            worst_probs = get_true_transition_probabilities(minPrs, maxPrs, collect(1:1:length(region_idx)))

            worst_exp_lb = sum(worst_probs .* sorted_actions[:, 3]') 

            sorted_actions2 = potential_actions[sortperm(potential_actions[:,4]), :] 
            region_idx = sorted_actions2[:,1]
            minPrs = sorted_actions2[:,5]
            maxPrs = sorted_actions2[:,6]
            worst_probs2 = get_true_transition_probabilities(minPrs, maxPrs, collect(1:1:length(region_idx))) 
            # worst_exp_ub = sum(worst_probs2 .* sorted_actions[:, 4]')  
            worst_exp_ub = sum(sorted_actions[:,4]')/length(sorted_actions[:,4]')

            if ocp.reward_type[1] == "worst-exp-lb" 
                total_reward = worst_exp_lb
            elseif ocp.reward_type[1] == "worst-exp-sum"
                sorted_actions2 = potential_actions[sortperm(potential_actions[:,4]), :] 
                region_idx = sorted_actions2[:,1]
                minPrs = sorted_actions2[:,5]
                maxPrs = sorted_actions2[:,6]
                worst_probs2 = get_true_transition_probabilities(minPrs, maxPrs, collect(1:1:length(region_idx)))
                total_reward = 0.5*(sum(worst_probs .* sorted_actions[:, 3]') + sum(worst_probs2 .* sorted_actions2[:, 4]'))
            elseif ocp.reward_type[1] == "worst-exp-lb-ub"
                # TOdo: Move all logic for this reward here
                total_reward = worst_exp_lb
            else
                @error "Unknown primary reward type!!!"
            end
          
            ###### Defining secondary reward just on DFA can result in weird behavior!!
            if ocp.reward_type[2] == "expected-dfa"
                x_new = [μ["x1"], μ["x2"]]
                x_new = reshape(x_new, 2, 1)
                secondary_reward = propogate_pimdp_trajectory_test(ocp.pimdp, ocp.dfa, ocp.extents, x_new)
            elseif ocp.reward_type[2] == "worst-dfa"
                secondary_reward = calculate_worst_dfa_reward(ocp.pimdp, ocp.dfa, imdp_extent_ids)
            elseif ocp.reward_type[2] == "best-dfa"
                secondary_reward = calculate_best_dfa_reward(ocp.pimdp, ocp.dfa, imdp_extent_ids)
            elseif ocp.reward_type[2] == "exp-probs"
                secondary_reward = (sum(minPrs) + sum(maxPrs))/length(maxPrs)
            elseif ocp.reward_type[2] == "combo"
                dfa_metric = calculate_worst_dfa_reward(ocp.pimdp, ocp.dfa, imdp_extent_ids)
                if dfa_metric > -Inf
                    secondary_reward = sum(minPrs) + sum(maxPrs)
                else
                    secondary_reward = dfa_metric
                end
            elseif ocp.reward_type[2] == "combo2"
                x_new = [μ["x1"], μ["x2"]]
                x_new = reshape(x_new, 2, 1)
                dfa_metric = propogate_pimdp_trajectory_test(ocp.pimdp, ocp.dfa, ocp.extents, x_new)
                if dfa_metric > -Inf
                    if isnothing(pimdp_reward_dict)
                        secondary_reward = calculate_pimdp_metric(ocp, mode)
                    else
                        secondary_reward = -pimdp_reward_dict[q_c]
                    end
                else
                    secondary_reward = dfa_metric
                end
            elseif ocp.reward_type[2] == "worst-exp-ub"
                sorted_actions = potential_actions[sortperm(potential_actions[:,4]), :] 
                region_idx = sorted_actions[:,1]
                minPrs = sorted_actions[:,5]
                maxPrs = sorted_actions[:,6]
                worst_probs = get_true_transition_probabilities(minPrs, maxPrs, collect(1:1:length(region_idx)))
                secondary_reward = sum(worst_probs .* sorted_actions[:, 4]') 
            elseif ocp.reward_type[2] == "worst-exp-ub-same"
                sorted_actions = potential_actions[sortperm(potential_actions[:,4]), :] 
                secondary_reward = sum(worst_probs .* sorted_actions[:, 4]') 
            elseif ocp.reward_type[2] == "none"
                secondary_reward = -Inf 
            elseif ocp.reward_type[2] == "exp-pimdp-metric1"
                # x_new = [μ["x1"], μ["x2"]]
                # x_new = reshape(x_new, 2, 1)
                # pimdp_state_ex = propogate_pimdp_trajectory_test(ocp.pimdp, ocp.dfa, ocp.extents, x_new)

                secondary_reward = 0
                for ps in nt_pimdp_states
                    val = pimdp_reward_dict[ps] 
                    if val < Inf
                        secondary_reward -= pimdp_reward_dict[ps] 
                    end
                end 
                secondary_reward *= 1/length(nt_pimdp_states)

                # if isnothing(pimdp_reward_dict)
                #     secondary_reward = calculate_pimdp_metric(ocp, mode)
                # else
                #     secondary_reward = -pimdp_reward_dict[pimdp_state_ex]
                # end
            elseif ocp.reward_type[2] == "worst-dfa-exp-pimdp-metric1"
                dfa_metric = calculate_worst_dfa_reward(ocp.pimdp, ocp.dfa, imdp_extent_ids)
                if dfa_metric > -Inf
                    # x_new = [μ["x1"], μ["x2"]]
                    # x_new = reshape(x_new, 2, 1)
                    # pimdp_state_ex = propogate_pimdp_trajectory_test(ocp.pimdp, ocp.dfa, ocp.extents, x_new)
                    secondary_reward = 0
                    for ps in nt_pimdp_states
                        val = pimdp_reward_dict[ps] 
                        if val < Inf
                            secondary_reward -= pimdp_reward_dict[ps] 
                        end
                    end 
                    secondary_reward *= 1/length(nt_pimdp_states)
                    # mean([-pimdp_reward_dict[ps] for ps in nt_pimdp_states])
                    # if isnothing(pimdp_reward_dict)
                    #     secondary_reward = calculate_pimdp_metric(ocp, mode)
                    # else
                    #     secondary_reward = -pimdp_reward_dict[pimdp_state_ex]
                    # end
                else
                    secondary_reward = dfa_metric
                end
            else
                @error "Unknown sec reward type!!!"
            end
            
            if total_reward > best_reward[1]
                best_mode = mode
                best_reward[1] = total_reward
                best_reward[2] = secondary_reward 
                best_gp_info_dict = gp_info_dict
                best_ub = worst_exp_ub
            elseif total_reward == best_reward[1]
                if worst_exp_ub > best_ub && ocp.reward_type[1] == "worst-exp-lb-ub" 
                    best_mode = mode
                    best_ub = worst_exp_ub 
                    best_reward[2] = secondary_reward 
                    best_gp_info_dict = gp_info_dict 
                elseif secondary_reward > best_reward[2]
                    best_reward[2] = secondary_reward
                    best_mode = mode  
                    best_gp_info_dict = gp_info_dict 
                elseif secondary_reward == best_reward[2]
                    # This chooses a random action among a final tie
                    if rand(1)[1] > 0.5 
                        best_reward[2] = secondary_reward
                        best_mode = mode  
                        best_gp_info_dict = gp_info_dict 
                    end
                end
            end
        end
    end
    return best_mode, best_reward, best_gp_info_dict
end

"""
Get pimdp states with non-trivial transition bounds
"""
function get_nontrivial_next_pimdp_states(pimdp, dfa, pimdp_state, mode, foo)
    P_row_idx = map_pimdp_state_to_Pmat_row(pimdp_state, length(dfa.states), length(pimdp.actions), mode)
    # @info pimdp_state
    prow_max = pimdp.Pmax[P_row_idx, :]

    # Find all states with interval bounds with non-zero UB
    regions_with_nonzero_ub = findall(x->x>0., prow_max) 
    return regions_with_nonzero_ub
end

"""
Get pimdp states with non-trivial transition bounds for all actions
"""
function get_nontrivial_next_pimdp_states(pimdp, dfa, pimdp_state, num_actions)
    pimdp_states = []
    for mode=1:num_actions
        pimdp_states = pimdp_states ∪ get_nontrivial_next_pimdp_states(pimdp, dfa, pimdp_state, mode, "foo") 
    end
    # Get unique states
    unique!(pimdp_states)
    return pimdp_states
end

"""
PIMDP Distance Metric
"""
function calculate_pimdp_metric(pimdp, pimdp_state, dfa, num_actions; max_depth=Inf, helper_dict=nothing)
    # Get the current pimdp state
    # pimdp_state = pimdp.state_history[end]

    # Get the labels that will induce a transition
    target_labels = get_next_dfa_labels(dfa, pimdp_state[2])
    sink_labels = get_dfa_sink_labels(dfa)

    distance = 0
    # From the current PIMDP state, get the next possible states
    pimdp_states = get_nontrivial_next_pimdp_states(pimdp, dfa, pimdp_state, 1, "foo")

    pimdp_states_tried = []
    searching_flag = true 

    best_dist = Inf
    while !isempty(pimdp_states) && searching_flag && distance < max_depth
        distance += 1
        next_states = []
        
        for state in pimdp_states
            imdp_state_idx = Int(pimdp_col_to_imdp_state(state, length(dfa.states)))
            pimdp_statep =  (imdp_state_idx, pimdp_state[2])
            if pimdp_statep ∈ keys(helper_dict)
                searching_flag = false
                subdistance = helper_dict[pimdp_statep] + distance
                if subdistance < best_dist
                    best_dist = subdistance
                end 
            elseif pimdp.imdp_label_dict[imdp_state_idx] ∈ target_labels
                searching_flag = false
                if distance < best_dist
                    best_dist = distance
                end
            end

            push!(pimdp_states_tried, state) 
            if searching_flag
                next_states = next_states ∪ get_nontrivial_next_pimdp_states(pimdp, dfa, pimdp_statep, num_actions)
            end
        end

        pimdp_states = setdiff(next_states, pimdp_states_tried)
        @assert isempty(pimdp_states ∩ pimdp_states_tried)
    end
        # For each next possible state, check the imdp label
    return best_dist
end

"""
Calculate Min Distance to Accepting state
"""
function calculate_pimdp_metric_all_offline(pimdp, dfa, num_actions; max_depth=Inf)
    # Get the current pimdp state
    pimdp_dist_dict = Dict()

    # Get the labels that will induce a transition
    states_finished = []
    sink_labels = get_dfa_sink_labels(dfa)
    for imdp_idx in keys(pimdp.imdp_label_dict), dfa_state=1:length(dfa.states)
    # for imdp_idx in [1], dfa_state=1:length(dfa.states)
        pimdp_state = (imdp_idx, dfa_state)
        imdp_label = pimdp.imdp_label_dict[imdp_idx]

        if pimdp_state[2] == dfa.sink_state || imdp_label ∈ sink_labels
            pimdp_dist_dict[pimdp_state] = Inf 
        elseif pimdp_state[2] == dfa.accepting_state || possible_accept_transition(dfa, dfa_state, imdp_label) 
            pimdp_dist_dict[pimdp_state] = 0
        elseif pimdp_state ∈ keys(pimdp_dist_dict)
            continue
        else
            val = calculate_pimdp_metric(pimdp, pimdp_state, dfa, num_actions, helper_dict=pimdp_dist_dict)
            pimdp_dist_dict[pimdp_state] = val
        end
        #     # From the current PIMDP state, get the next possible states
        #     pimdp_states = get_nontrivial_next_pimdp_states(pimdp, dfa, pimdp_state, num_actions)
        #     @info pimdp_states
        #     pimdp_states_tried = []
        #     searching_flag = true 
        #     distance = 0
        #     while !isempty(pimdp_states) && searching_flag && distance < max_depth
        #         distance += 1
        #         @info distance
        #         next_states = []
        #         for state in pimdp_states

        #             imdp_idxp = pimdp_col_to_imdp_state(state, length(dfa.states))
        #             imdp_labelp = pimdp.imdp_label_dict[imdp_idxp]
        #             dfa_statep = state % length(dfa.states) 
        #             pimdp_statep = (Int(imdp_idxp), dfa_statep)

        #             if (pimdp_statep[2] == 0 && imdp_labelp ∉ sink_labels) || possible_accept_transition(dfa, dfa_statep, imdp_labelp) 
        #                 @info "STATE FOUND --------"
        #                 @info pimdp_statep
        #                 searching_flag = false
        #                 pimdp_dist_dict[pimdp_state] = distance 
        #             elseif pimdp_statep in keys(pimdp_dist_dict) 
        #                 searching_flag = false
        #                 pimdp_dist_dict[pimdp_state] = pimdp_dist_dict[pimdp_statep] + distance
        #             end
        #             push!(pimdp_states_tried, state) 
        #             next_states = next_states ∪ get_nontrivial_next_pimdp_states(pimdp, dfa, pimdp_statep, num_actions)
        #             @info next_states
        #         end
        #         pimdp_states = setdiff(next_states, pimdp_states_tried)
        #         # Check for states already procced here
        #         @assert isempty(pimdp_states ∩ pimdp_states_tried)
        #     end
        # end
        @assert pimdp_state ∈ keys(pimdp_dist_dict) 
    end
    return pimdp_dist_dict 
end


function get_distance_to_dfa_label(ocp, mode)
    # TODO: Finish this
   # Get the current pimdp state
   pimdp_state = ocp.pimdp.state_history[end]

   # Get the labels that will induce a transition
   target_labels = get_next_dfa_labels(ocp.dfa, pimdp_state[2])

   # From the current PIMDP state, get the next possible states
   pimdp_states = get_nontrivial_next_pimdp_states(ocp, pimdp_state, mode)

   distance = Inf
   for state in pimdp_states
        imdp_state_idx = Int(pimdp_col_to_imdp_state(state, length(ocp.dfa.states))) 

   end
end

"""
Get the policy from a synthesis result saved as a matrix.
"""
function get_policy_from_offline_synthesis_mat(filename)
    policy = Dict()
    synth_dict = matread(filename)

    policy = synth_dict["policy"]
    minPr = synth_dict["indVmin"]
    maxPr = synth_dict["indVmax"]

    for (i, action) in enumerate(synth_dict["policy"])
        policy[i] = Int(action)
    end
    return policy, minPr, maxPr
end

