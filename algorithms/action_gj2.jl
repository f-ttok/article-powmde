using Arpack
using FastGaussQuadrature
using LinearAlgebra
using SparseArrays
using SuiteSparse


function powmvec_gj2_spd(A, b, α; ϵ=1e-6, m_max=1000)
    λmax = eigs(A, nev=1, tol=1e-3, which=:LM)[1][1]
    λmin = eigs(A, nev=1, tol=1e-3, which=:SM)[1][1]
    c = 1 / sqrt(λmax*λmin)

    A_tilde = c * A
    f(v) = ((1-v)*I + (1+v)*A_tilde) \ b
    λmax_tilde = c * λmax
    ϵ_tilde = c^α * ϵ


    m_min = 2
    m = (m_max + m_min) ÷ 2
    abserr = Inf
    exact = λmax_tilde^α
    while m_max - m_min > 1
        m = (m_max + m_min) ÷ 2
        approx = pow_gj2(λmax_tilde, α, m)
        abserr = abs(exact - approx)
        if abserr > ϵ_tilde
            m_min = m
        else
            m_max = m
        end
    end
    m = m_max

    v, w = gaussjacobi(m, α-1, -α)
    t = zero(b)
    for k = 1:m
        t .+= w[k] * f(v[k])
    end
    t .= c^(-α) * 2sin(α*π) / π * (A_tilde * t)
    return t, m
end


function pow_gj2(λ, α, m)
    v, w = gaussjacobi(m, α-1, -α)
    approx = sum(@. w / ((1-v) + (1+v)*λ))
    approx *= 2sin(α*π) / π * λ
    return approx
end


function powmvec_gj2(A, b, α; ϵ=1e-6, m0=8, eval_max=1000)
    σ_max = svds(A, nsv=1, tol=1e-3)[1].S[1]
    σ_min = sqrt(real(eigs(A'*A, nev=1, tol=1e-3, which=:SM)[1][1]))
    c = 1 / sqrt(σ_max*σ_min)

    A_tilde = c * A
    f(v) = ((1-v)*I + (1+v)*A_tilde) \ b
    ϵ_tilde = c^α * ϵ

    m = m0
    nof_eval = m0
    t_old = zero(b)
    v, w = gaussjacobi(m, α-1, -α)
    for k = 1:m
        t_old .+= w[k] * f(v[k])
    end
    approx_old = 2sin(α*π) / π * (A_tilde * t_old)

    t_new = zero(b)
    approx_new = zero(b)
    err_est = Inf
    while err_est > ϵ_tilde && nof_eval < eval_max
        m = 2m
        nof_eval += m
        v, w = gaussjacobi(m, α-1, -α)
        t_new .= 0
        for k = 1:m
            t_new .+= w[k] * f(v[k])
        end
        approx_new = 2sin(α*π) / π * (A_tilde * t_new)
        err_est = norm(approx_new - approx_old)

        t_old .= t_new
        approx_old .= approx_new
    end
    if err_est > ϵ_tilde
        println("warning: the algorithm did not converged (err_est = $(err_est))")
    end

    t_new .=  c^(-α) * 2sin(α*π) / π * (A_tilde * t_new)
    return t_new, nof_eval
end