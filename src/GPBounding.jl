module GPBounding

using GaussianProcesses
using LinearAlgebra

include("squared_exponential.jl")

export compute_μ_bounds_bnb, compute_σ_ub_bounds, compute_σ_ub_bounds_auto
#TODO: Eventually move this to a different file. 
#TODO: Compute lower bounds of sigma someday

function compute_μ_bounds_bnb(gp, x_L, x_U; max_iterations=100, bound_epsilon=1e-2, max_flag=false)
    # By default, it calculates bounds on the minimum. 
    # Calculate this vector outside of the loop to save computation.
    theta_vec_train_squared = zeros(gp.nobs);
    theta_vec = ones(gp.dim) * 1 ./ (2*gp.kernel.ℓ2)
    for i = 1:gp.nobs
        theta_vec_train_squared[i] = transpose(theta_vec) * (gp.x[:, i].^2)
    end   
    
    x_best, lbest, ubest = compute_μ_lower_bound(gp, x_L, x_U, theta_vec_train_squared, upper_flag=max_flag)
    if max_flag
        temp = lbest
        lbest = -ubest
        ubest = -temp
    end
    
    candidates = [[(x_L, x_U), lbest, ubest]]
    iterations = 0
    lb_list = []
   
    while !isempty(candidates) && iterations < max_iterations
        new_candidates = []
        for candidate in candidates
            
            extent = candidate[1]
            lb_can = candidate[2]
            ub_can = candidate[3]
            bound_pairs = split_region(extent[1], extent[2])
            
            for pair in bound_pairs
                x_lb1, lb1, ub1 = compute_μ_lower_bound(gp, pair[1], pair[2], theta_vec_train_squared, upper_flag=max_flag)

                if max_flag
                    temp = lb1
                    lb1 = -ub1
                    ub1 = -temp
                end
                
                if ub1 <= ubest
                    ubest = ub1
                    lbest = lb1
                    x_best = x_lb1
                    push!(new_candidates, hcat([pair, lb1, ub1])) 
                elseif lb1 < ubest   
                    push!(new_candidates, hcat([pair, lb1, ub1]))    
                end
                
            end
            
        end

        if norm(ubest - lbest) < bound_epsilon
            if max_flag
                temp = lbest
                lbest = -ubest
                ubest = -temp
            end
            return x_best, lbest, ubest
        end
        candidates = new_candidates
        iterations += 1
    end
    if max_flag
        temp = lbest
        lbest = -ubest
        ubest = -temp
    end
    return x_best, lbest, ubest  
end

function compute_σ_ub_bounds(gp, K_inv, x_L, x_U; max_iterations=10, bound_epsilon=1e-4)
    
    lbest = -Inf
    ubest = 1.
    x_best = nothing
    
    candidates = [[(x_L, x_U), lbest, ubest]]
    iterations = 0
   
    while !isempty(candidates) && iterations < max_iterations
        new_candidates = []
        for candidate in candidates
            
            extent = candidate[1]
            lb_can = candidate[2]
            ub_can = candidate[3]

            if ub_can < lbest
                continue
            end
            
            round_lbest = lbest
            round_ubest = ubest
            
            bound_pairs = split_region(extent[1], extent[2])

            for pair in bound_pairs
                x_ub1, lb1, ub1 = compute_σ_upper_bound(gp, pair[1], pair[2], K_inv)
                
                if lb1 > lbest
                    lbest = lb1
                    ubest = ub1
                    x_best = x_ub1
                    push!(new_candidates, hcat([pair, lb1, ub1]))
                elseif ub1 < ubest && ub1 > lbest
                    ubest = ub1
                    push!(new_candidates, hcat([pair, lb1, ub1]))
                elseif lb1 > lb_can && ub1 < ub_can && ub1 > lbest
                    push!(new_candidates, hcat([pair, lb1, ub1]))
                end
                
                if norm(ubest - lbest) < bound_epsilon
                    @debug ubest, lbest
                    return x_best, lbest, ubest
                end
            end

        end
        
        if (length(new_candidates) == 0 && norm(ubest - lbest) > bound_epsilon)
            delta_x = 10. ^(-iterations)
            if delta_x < 1e-6
                break
            end
            x_L_n = [x_best[i] - delta_x for i=1:length(x_best)]
            x_U_n = [x_best[i] + delta_x for i=1:length(x_best)]
            new_candidates = [[(x_L_n, x_U_n), lbest, ubest]]
        end
        candidates = new_candidates
        iterations += 1
    end
    
    while norm(ubest - lbest) > bound_epsilon
        delta_x = 10. ^(-iterations)
        if delta_x < 1e-6
            break
        end
        x_L_n = [x_best[i] - delta_x for i=1:length(x_best)]
        x_U_n = [x_best[i] + delta_x for i=1:length(x_best)]
        x_best, lbest, ubest = compute_σ_upper_bound(gp, x_L_n, x_U_n, K_inv)
        iterations += 1
    end
    
    return x_best, lbest, ubest
    
end

end # module
