using Arpack
using LinearAlgebra
using SparseArrays
using SuiteSparse


function powmvec_de_spd(A, b, α; ϵ=1e-6, m_max=1000)
    Δ = 1e-3
    λ_max = eigs(A, nev=1, tol=Δ, which=:LM)[1][1]
    λ_min = eigs(A, nev=1, tol=Δ, which=:SM)[1][1]
    c = 1 / sqrt(λ_max*λ_min)

    A_tilde = c * A
    f(x) = exp(α*π*sinh(x)/2) * cosh(x) * ((exp(π*sinh(x)/2)*I+A_tilde) \ b)
    λ_max_tilde = c * λ_max
    ϵ_tilde = c^α * ϵ

    norm_A, norm_A_inv = λ_max_tilde, λ_max_tilde
    l, r = get_interval(norm_A, norm_A_inv, α, ϵ_tilde/(1 + 1/(1-Δ)))
    m_min = 2
    m = (m_max + m_min) ÷ 2
    abserr = Inf
    exact = λ_max_tilde^α

    while m_max - m_min > 1
        m = (m_max + m_min) ÷ 2
        approx = pow_de(λ_max_tilde, α, m, l, r)
        abserr = abs(exact - approx)
        if abserr > ϵ_tilde
            m_min = m
        else
            m_max = m
        end
    end
    m = m_max

    x = LinRange(l, r, m)
    h = (r-l)/(m-1)
    w = fill(h, m)
    w[1], w[m] = h/2, h/2

    t1, t2 = zero(b), zero(b)
    for k = 1:m÷2
        t1 .+= w[k] * f(x[k])
    end
    for k = m:-1:m÷2+1
        t2 .+= w[k] * f(x[k])
    end
    t1 .= c^(-α) * sin(α*π) / 2 * (A_tilde * (t1 + t2))
    return t1, m
end


function get_interval(norm_A, norm_A_inv, α, ϵ)
    a1 = ϵ*π*α*(1+α) / (4*sin(α*π)*(1+2*α))
    a2 = (2*norm_A_inv)^(-α)
    a = min(a1, a2)
    b1 = (ϵ*π*(1-α)*(2-α))^(α/(α-1)) / (4*sin(α*π)*(3-2α)*norm_A)^(α/(α-1))
    b2 = (2*norm_A)^α
    b = max(b1, b2)
    l = asinh(2*log(a)/(α*π))
    r = asinh(2*log(b)/(α*π))
    return l, r
end


function pow_de(λ, α, m, l, r)
    f(x) = exp(α*π*sinh(x)/2) * cosh(x) / (exp(π*sinh(x)/2)+λ)
    x = LinRange(l, r, m)
    h = (r-l)/(m-1)
    w = fill(h, m)
    w[1], w[m] = h/2, h/2

    t1, t2 = 0.0, 0.0
    for k = 1:m÷2
        t1 += w[k] * f(x[k])
    end
    for k = m:-1:m÷2+1
        t2 += w[k] * f(x[k])
    end
    return sin(α*π)/2 * λ * (t1 + t2)
end



function powmvec_de(A, b, α; ϵ=1e-6, m0=8, m_max=1000)
    Δ = 1e-3
    σ_max = svds(A, nsv=1, tol=Δ)[1].S[1]
    σ_min = sqrt(real(eigs(A'*A, nev=1, tol=Δ, which=:SM)[1][1]))

    c = 1 / sqrt(σ_max*σ_min)

    A_tilde = c * A
    f(x) = exp(α*π*sinh(x)/2) * cosh(x) * ((exp(π*sinh(x)/2)*I+A_tilde) \ b)
    ϵ_tilde = c^α * ϵ
    σ_max_tilde = c * σ_max
    σ_min_tilde = c * σ_min

    norm_A, norm_A_inv = σ_max_tilde, 1/σ_min_tilde
    l, r = get_interval(norm_A, norm_A_inv, α, ϵ_tilde/(1 + 1/(1-Δ)))

    m = m0
    x = LinRange(l, r, m)
    h = (r-l) / (m-1)
    w = fill(h, m)
    w[1], w[m] = h/2, h/2
    t1, t2 = zero(b), zero(b)
    for k = 1:m÷2
        t1 .+= w[k] * f(x[k])
    end
    for k = m:-1:m÷2+1
        t2 .+= w[k] * f(x[k])
    end
    t_old = t1 + t2
    approx_old = sin(α*π) / 2 * (A_tilde * t_old)

    t_new = zero(t_old)
    err_est = Inf
    while err_est > ϵ_tilde && m < m_max
        h = h/2
        t1 .= 0
        t2 .= 0
        for k = 1:(m-1)÷2
            t1 .+= h * f(l + (2k-1)*h)
        end
        for k = (m-1):-1:((m-1)÷2+1)
            t2 .+= h * f(l + (2k-1)*h)
        end
        t_new .= t_old/2 + (t1 + t2)
        approx_new = sin(α*π) / 2 * (A_tilde * t_new)
        err_est = norm(approx_new - approx_old)
        shoud_stop =
        m = 2m - 1
        t_old .= t_new
        approx_old .= approx_new
    end
    if err_est > ϵ_tilde
        println("warning: the algorithm did not converged (err_est = $(err_est))")
    end

    t_new .=  c^(-α) * sin(α*π) / 2 * (A_tilde * t_new)
    return t_new, m
end
