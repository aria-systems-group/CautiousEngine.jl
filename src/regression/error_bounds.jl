"Calculate the value of the Munich bound in Lederer et. al.'s paper (add reference)."
function munich_bound(point, gp, prior_norm, delta; Lf=2., Xradius=10.)
    # Get the value of the standard deviation at the test point
    _, σ_sq= predict_y(gp, point)
    σ = sqrt(σ_sq[1])
    # Get the total number of points used for training
    N = length(gp.y)
    # Calculate a suitable grid discretization
    tau = 1e-6/(N^2)
    # Get the kernel parameters used in the bound
    sf = gp.kernel.σ2
    ls = gp.kernel.ℓ2 
    kernel_max = sf 
    # Calculate the Lipschitz constant of the kernel
    Lk = norm(sf*exp(-0.5)./sqrt(ls))
    # Calculate norms used in the bound
    InvNormY = norm(gp.alpha) 
    InvNorm = prior_norm 
    # Calculate the bound of the mean function Lipschitz constant
    Lv_bound = Lk*sqrt(N)*InvNormY
    # Calculate the bound on the covering number 
    M_bound = (1+Xradius/tau)^2 
    # Calculate the bound on omega
    omega_T_bound = sqrt(2*tau*Lk*(1+N*InvNorm*kernel_max)) 
    # Calculate Beta
    Beta_T = 2*log(M_bound/delta)
    # Calculate Gamma
    gamma_T = (Lv_bound + Lf)*tau +sqrt(Beta_T)*omega_T_bound
    # Put it all together to get the bound
    bound = sqrt(Beta_T)*σ + gamma_T
    return bound
end

" Calculate an upper bound on the probability parameter delta given an epsilon."
function munich_bound_prob(point, gp, prior_norm, epsilon; Lf=2., Xradius=10.)
    # Get the value of the standard deviation at the test point
    _, σ2 = predict_y(gp, point)
    σ = sqrt(σ2[1])
    # Get the total number of points used for training
    N = length(gp.y)
    # Calculate a suitable grid discretization
    tau = 1e-6/(N^2)
    # Get the kernel parameters used in the bound
    sf = gp.kernel.σ2
    ls = gp.kernel.ℓ2 
    kernel_max = sf 
    # Calculate the Lipschitz constant of the kernel
    Lk = norm(sf*exp(-0.5)./sqrt(ls))
    # Calculate norms used in the bound
    InvNormY = norm(gp.alpha) 
    InvNorm = prior_norm 
    # Calculate the bound of the mean function Lipschitz constant
    Lv_bound = Lk*sqrt(N)*InvNormY
    # Calculate the bound on the covering number 
    M_bound = (1+Xradius/tau)^2 
    # Calculate the bound on omega
    omega_T_bound = sqrt(2*tau*Lk*(1+N*InvNorm*kernel_max)) 
    delta_UB = M_bound*exp(-0.5*(epsilon-(Lv_bound+Lf)*tau)^2/(σ + omega_T_bound)^2)
    delta_UB = minimum([delta_UB, 1.])
    return delta_UB
end

# TODO: The values for the following two functions are heuristic. They are not guaranteeed to be good. 
"Calculate the delta-confidence interval using the RKHS-approach (add reference)"
function original_rkhs_bound(σ, δ, gp_info)
    B = gp_info.RKHS_bound
    γ_bd = gp_info.γ_bound
    T = gp_info.gp.nobs
    C = sqrt(2*B^2 + 300*γ_bd*log((T+1)/δ)^3)
    return σ*C 
end

"Calculate an upper bound on the probability parameter given a value of epsilon."
function original_rkhs_bound_prob(σ, ϵ, gp_info)
    γ_bd = gp_info.γ_bound
    B = gp_info.RKHS_bound
    T = gp_info.gp.nobs
    dbound = exp(-cbrt(((ϵ/σ)^2 - 2*B^2) / (300*γ_bd)) + log(T+1))
    return minimum([dbound, 1.])
end

"New RKHS Bound from Theorem 2 of Chowdhury et al"
function tight_rkhs_bound(σ, δ, gp_info)
    γ_bd = gp_info.γ_bound 
    B = gp_info.RKHS_bound
    R = gp_info.post_scale_factor*exp(gp_info.gp.logNoise.value)
    C = gp_info.post_scale_factor*(B + R*sqrt(2*(γ_bd + 1 + log(1. /δ))))
    return σ*C
end

"New RKHS Bound from Theorem 2 of Chowdhury et al"
function tight_rkhs_bound_prob(σ, ϵ, gp_info)
    γ_bd = gp_info.γ_bound
    B = gp_info.RKHS_bound
    R = gp_info.post_scale_factor*exp(gp_info.gp.logNoise.value)

    frac = ϵ/(gp_info.post_scale_factor*σ)
    if frac > B
        dbound = exp(-0.5*(1/R*(frac - B))^2 + γ_bd + 1.)
    else
        dbound = 1. 
    end
    return minimum([dbound, 1.])
end

