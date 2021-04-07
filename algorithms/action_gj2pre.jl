using Arpack
using FastGaussQuadrature
using GSL
using LinearAlgebra
using SparseArrays
using SuiteSparse


"""
    τ(λmax, λmin, α, m)

Selecting a parameter τ>0 for GJ2 [1, Eq. (37)].

[1] Aceto, L., & Novati, P. (2019). Rational approximations to fractional powers of self-adjoint positive operators. Numer. Math., 143(1), 1--16. https://doi.org/10.1007/s00211-019-01048-4

"""
function τ(λmax, λmin, α, m)
    μmax, μmin = 1/λmin, 1/λmax
    κ = μmax / μmin
    e = exp(1)
    m_bar = α / (2sqrt(2)) * sqrt(log(κ * e^2)) * κ^(1/4)
    
    if m < m_bar
        w = sf_lambert_W0(4*m^2*e/α^2)
        t = μmin*(α/2/m/e)^2 * exp(2*w)
        return t
    else
        tmp = α*sqrt(μmax)/8/m * log(κ)
        t = (-tmp + sqrt(tmp^2 + sqrt(μmax*μmin)))^2
        return t
    end
end


"""
    errest_gj2plus(λmax, λmin, α, m)

Error estimate.
"""
function errest_gj2pre(λmax, λmin, α, m)
    μmax, μmin = 1/λmin, 1/λmax
    κ = λmax / λmin
    e = exp(1)
    m_bar = α / (2sqrt(2)) * sqrt(log(κ * e^2)) * κ^(1/4)
    
    if m < m_bar
        est = 2*sin(α*π)*μmin^-α
        est *= (2*m*sqrt(e)/α)^(-4α)
        est *= (2*log(2m/α)+1)^(2α)
        return est
    else
        est = 2*sin(α*π)*(μmin*μmax)^(-α/2)
        est *= exp(-4m*κ^(-1/4))
        return est
    end
end


function powmvec_gj2pre_spd(A, b, α; ϵ=1e-6, m_max=1000)
    λmax = eigs(A, nev=1, tol=1e-3, which=:LM)[1][1]
    λmin = eigs(A, nev=1, tol=1e-3, which=:SM)[1][1]

    f(v, A_tilde) = ((1-v)*I + (1+v)*A_tilde) \ b

    m_min = 2
    m = (m_max + m_min) ÷ 2
    abserr = Inf
    while m_max - m_min > 1
        m = (m_max + m_min) ÷ 2
        abserr = errest_gj2pre(λmax, λmin, α, m)
        if abserr > ϵ
            m_min = m
        else
            m_max = m
        end
    end
    m = m_max

    v, w = gaussjacobi(m, α-1, -α)
    t = zero(b)
    c = τ(λmax, λmin, α, m)
    A_tilde = c * A
    for k = 1:m
        t .+= w[k] * f(v[k], A_tilde)
    end
    t .= c^(-α) * 2sin(α*π) / π * (A_tilde * t)
    return t, m
end