# Methods related to online control of dynamic systems

"""
Objects relevant to online control problem
"""
struct OnlineControlProblem
    initial_state
    control_modes
    horizon::Int
    imdp::IMDP
    pimdp::PIMDP
    dfa::DFA
    extents
    control_policy
    local_gps::Dict
end

"""
Setup the online control problem.
"""
function setup_online_control_problem(x0, modes, steps, imdp, pimdp, dfa, extents, policy, local_gps)
    reset_pimdp(x0, imdp, dfa, pimdp, extents, policy)
    return OnlineControlProblem(x0, modes, steps, imdp, pimdp, dfa, extents, policy, local_gps)
end

"""
One step of a control loop with an offline policy.
"""
function online_control_loop(ocp::OnlineControlProblem)
    x_c = ocp.pimdp.trajectory[end]
    mode = ocp.pimdp.policy[map_pimdp_state_to_index(ocp.pimdp, ocp.pimdp.state_history[end])]
    x_n = ocp.modes[Int(mode)](x_c)
    propogate_pimdp_trajectory(ocp.pimdp, ocp.dfa, ocp.extents, ocp.x_n)
end

"""
N-Step control loop with an offline policy.
"""
function online_control_loop(ocp::OnlineControlProblem, steps::Int)
    total_time = 0.
    for i=1:steps
        etime = @elapsed begin
            # check_pimdp_exit_conditions(ocp.pimdp) ? break : nothing
            x_c = ocp.pimdp.trajectory[end]
            q_c = ocp.pimdp.state_history[end]
            pimdp_state = (q_c[1]*length(ocp.dfa.states)) - length(ocp.dfa.states) + q_c[2]
            mode = ocp.control_policy[pimdp_state]
            # mode = ocp.pimdp.policy[map_pimdp_state_to_index(ocp.pimdp, ocp.pimdp.state_history[end])]
            x_n = ocp.control_modes[mode](x_c)
            propogate_pimdp_trajectory(ocp.pimdp, ocp.dfa, ocp.extents, x_n)
        end
        @info "Time in loop: $etime seconds"
        total_time += etime
    end
    @info "Total time for $steps steps: $total_time seconds"
end

"""
One step optimal control loop with an offline policy and GP update.
"""
function online_control_loop_optimal(ocp::OnlineControlProblem)
    etime = @elapsed begin
        # check_pimdp_exit_conditions(ocp.pimdp) ? break : nothing
        x_c = ocp.pimdp.trajectory[end]
        q_c = ocp.pimdp.state_history[end]

        pimdp_state = (q_c[1]*length(ocp.dfa.states)) - length(ocp.dfa.states) + q_c[2]
        mode = ocp.control_policy[pimdp_state]
        # mode = ocp.pimdp.policy[map_pimdp_state_to_index(ocp.pimdp, ocp.pimdp.state_history[end])]
        x_n = ocp.control_modes[mode](x_c)
        propogate_pimdp_trajectory(ocp.pimdp, ocp.dfa, ocp.extents, x_n)

        # GP Directory
        gp_dir = ocp.local_gps[mode]
        imdp_state = q_c[1]
        @info gp_dir, imdp_state
        gp_filename = filter(x -> occursin("-$imdp_state.bin", x), readdir(gp_dir))[1]
        gp_info = load_gp("$gp_dir/$gp_filename")
        gp_set = gp_info[:gp_set]
        for (i, gp_key) in enumerate(keys(gp_set))
            gp = gp_set[gp_key]
            y_meas = x_c[i] - x_n[i]
            update_gp(gp, (x_c, y_meas))
        end
    end
    @info "Total time: $etime seconds"
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
