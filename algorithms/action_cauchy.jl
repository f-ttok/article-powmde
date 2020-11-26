using SparseArrays
using SuiteSparse
using Elliptic

function powmvec_cauchy_spd(A, b, α; ϵ=1e-6, m_max=1000)
    # Compute the extreme eigenvalues
    λmax = eigs(A, nev=1, tol=1e-3, which=:LM)[1][1]
    λmin = eigs(A, nev=1, tol=1e-3, which=:SM)[1][1]
    c = 1 / sqrt(λmax*λmin)

    A_tilde = c * A
    λmax_tilde = c * λmax
    λmin_tilde = c * λmin
    ϵ_tilde = c^α * ϵ

    # Select the number of abscissas
    m_min = 2
    m = (m_max + m_min) ÷ 2
    abserr = Inf
    exact = [λmax_tilde^α, λmax_tilde^(-α)]

    while m_max - m_min > 1
        m = (m_max + m_min) ÷ 2
        approx = pow_cauchy(λmax_tilde, α, m)
        abserr = norm(exact - approx, Inf)
        if abserr > ϵ_tilde
            m_min = m
        else
            m_max = m
        end
    end
    m = m_max

    κ = λmax/λmin
    k = (κ^(1/4) - 1) / (κ^(1/4) + 1)
    K = Elliptic.K(k^2)
    Kp = Elliptic.K(1 - k^2)

    t = collect(1:m) .- 0.5
    t = @. im/2*Kp - K + t*2*K/m
    u, cn, dn = zero(t), zero(t), zero(t)
    for i = 1:m
        u[i], cn[i], dn[i] = ellipj(t[i], k^2)
    end

    w = @. (λmax_tilde*λmin_tilde)^(1/4) * (1/k+u) / (1/k-u)
    dzdt = @. cn * dn / (1/k-u)^2

    f = zeros(ComplexF64, size(b))
    for j = 1:m
        f .+= w[j]^(2α-1) * dzdt[j] * ((w[j]^2*I - A_tilde) \ b)
    end
    f .= c^(-α) * -8*K*(λmax_tilde*λmin_tilde)^(1/4) * (A_tilde * imag.(f)) / (k*π*m)
    return f, m
end


function pow_cauchy(λmax, α, m)
    λmin = 1 / λmax
    A = Diagonal([λmax, λmin])
    κ = λmax / λmin
    k = (κ^(1/4) - 1) / (κ^(1/4) + 1)
    K = Elliptic.K(k^2)
    Kp = Elliptic.K(1 - k^2)

    t = collect(1:m) .- 0.5
    t = @. im/2*Kp - K + t*2*K/m
    u, cn, dn = zero(t), zero(t), zero(t)
    for i = 1:m
        u[i], cn[i], dn[i] = ellipj(t[i], k^2)
    end

    w = @. (λmax*λmin)^(1/4) * (1/k+u) / (1/k-u)
    dzdt = @. cn * dn / (1/k-u)^2

    S = zeros(ComplexF64, 2, 2)
    for j = 1:m
        S += w[j]^(2α-1) * inv(w[j]^2*I - A) * dzdt[j]
    end
    S = -8*K*(λmax*λmin)^(1/4) * imag.(S) * A / (k*π*m)
    return real.(diag(S))
end


ellipj(u::T, m) where {T<:Real} = Elliptic.ellipj(u, m)


function ellipj(z::T, m) where {T<:Complex}
    x = real(z)
    y = imag(z)

    s, c, d = ellipj(x, m)
    s1, c1, d1 = ellipj(y, 1-m)

    denom = c1^2 + m*s^2*s1^2
    sn = (s*d1 + im*c*d*s1*c1) / denom
    cn = (c*c1 - im*s*d*s1*d1) / denom
    dn = (d*c1*d1 - im*m*s*c*s1) / denom

    return sn, cn, dn
end