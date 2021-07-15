module GPBounding

using GaussianProcesses
using LinearAlgebra
using Random
using Distributions

include("squared_exponential.jl")

export compute_μ_bounds_bnb, compute_σ_ub_bounds, compute_σ_ub_bounds_auto, compute_σ_ub_bounds_approx, compute_σ_ub_bounds_from_gp

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

function compute_σ_ub_bounds_approx(gp, x_L, x_U)
    N = 100
    mt = MersenneTwister(11)
    σ2_best = -Inf
    # Get N samples uniformly dist.
    x_samp = zeros(2,1)
    # @info x_samp
    for i=1:N
        # x_samp = vcat([rand(mt, Uniform(x_L[1], x_U[1]), 1, 1)[1], rand(mt, Uniform(x_L[2], x_U[1]), 1, 1)[1]])
        x_samp[1] = rand(mt, Uniform(x_L[1], x_U[1]), 1, 1)[1]
        x_samp[2] = rand(mt, Uniform(x_L[2], x_U[2]), 1, 1)[1]
        # @info x_samp
        _, σ2 = predict_f(gp, x_samp)
        if σ2[1] > σ2_best
            σ2_best = σ2[1] 
        end
    end
    return 0., 0., sqrt(σ2_best[1])
end

function create_x_matrix(xL, xU, N)
    dim = length(xL)
    xds = [(xU[i]-xL[i])/N for i=1:dim]
    # Assume 2D for now
    xd_prod = collect(Iterators.product(xL[1]:xds[1]:xU[1], xL[2]:xds[2]:xU[2]))
    return reshape(collect(Iterators.flatten(xd_prod)), dim, length(xd_prod))
end

function prepare_σ_gp(gp, x_L, x_U, ub; N=10, ll=-2.0)
    xn = create_x_matrix(x_L, x_U, N)
    yn = sqrt.(predict_f(gp, xn)[2])
    
    meanfcn = MeanConst(ub)
    kernel = SE(ll, 0.0)
    
    gp_σ = GP(xn, yn, meanfcn, kernel, 0.0)
    return gp_σ
end

function compute_σ_ub_bounds_from_gp(gp, x_L, x_U; ub=1.0)
    σ_gp = prepare_σ_gp(gp, x_L, x_U, ub)
    res_test = GPBounding.compute_μ_bounds_bnb(σ_gp, x_L, x_U, max_flag=true, max_iterations=4)
    return res_test[1], 0., res_test[3][1]+ub
end

end # module
