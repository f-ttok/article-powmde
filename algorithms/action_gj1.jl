using Arpack
using FastGaussQuadrature
using LinearAlgebra
using SparseArrays
using SuiteSparse


function powmvec_gj1_spd(A, b, α; ϵ=1e-6, m_max=1000)
    λ_max = eigs(A, nev=1, tol=1e-3, which=:LM)[1][1]
    λ_min = eigs(A, nev=1, tol=1e-3, which=:SM)[1][1]
    c = 1 / sqrt(λ_max*λ_min)

    A_tilde = c * A
    f(u) = ((1+u)^(1/α)*I + (1-u)^(1/α)*A_tilde) \ b
    λ_max_tilde = c * λ_max
    ϵ_tilde = c^α * ϵ

    m_min = 2
    m = (m_max + m_min) ÷ 2
    abserr = Inf
    exact = λ_max_tilde^α
    while m_max - m_min > 1
        m = (m_max + m_min) ÷ 2
        approx = pow_gj1(λ_max_tilde, α, m)
        abserr = abs(exact - approx)
        if abserr > ϵ_tilde
            m_min = m
        else
            m_max = m
        end
    end
    m = m_max

    u, w = gaussjacobi(m, 1/α-2, 0)
    t = zero(b)
    for k = 1:m
        t .+= w[k] * f(u[k])
    end
    t .= c^(-α) * 2sin(α*π) / (α*π) * (A_tilde * t)
    return t, m
end


function pow_gj1(λ, α, m)
    u, w = gaussjacobi(m, 1/α-2, 0)
    approx = sum(@. w / ((1+u)^(1/α) + (1-u)^(1/α)*λ))
    approx *= 2sin(α*π) / (α*π) * λ
    return approx
end


function powmvec_gj1(A, b, α; ϵ=1e-6, m0=8, eval_max=1000)
    σ_max = svds(A, nsv=1, tol=1e-3)[1].S[1]
    σ_min = sqrt(real(eigs(A'*A, nev=1, tol=1e-3, which=:SM)[1][1]))
    c = 1 / sqrt(σ_max*σ_min)

    A_tilde = c * A
    f(u) = ((1+u)^(1/α)*I + (1-u)^(1/α)*A_tilde) \ b
    ϵ_tilde = c^α * ϵ

    m = m0
    nof_eval = m0
    t_old = zero(b)
    u, w = gaussjacobi(m, 1/α-2, 0)
    for k = 1:m
        t_old .+= w[k] * f(u[k])
    end
    approx_old = 2sin(α*π) / (α*π) * (A_tilde * t_old)

    t_new = zero(b)
    approx_new = zero(b)
    err_est = Inf
    while err_est > ϵ_tilde && nof_eval < eval_max
        m = 2m
        nof_eval += m
        u, w = gaussjacobi(m, 1/α-2, 0)
        t_new .= 0
        for k = 1:m
            t_new .+= w[k] * f(u[k])
        end
        approx_new = 2sin(α*π) / (α*π) * (A_tilde * t_new)

        err_est = norm(approx_new - approx_old)
        t_old .= t_new
        approx_old .= approx_new
    end
    if err_est > ϵ_tilde
        println("warning: the algorithm did not converged (err_est = $(err_est))")
    end

    t_new .= c^(-α) * 2sin(α*π) / (α*π) * (A_tilde * t_new)
    return t_new, nof_eval
end