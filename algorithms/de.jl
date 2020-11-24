using LinearAlgebra


"""
    get_interval(norm_A, norm_A_inv, α, ϵ)

Computing a finite interval satisfying (3.1)
"""
function get_interval(norm_A, norm_A_inv, α, ϵ)
    a1 = ϵ*π*α*(1+α) / (4*sin(α*π)*(1+2*α))
    a2 = (2*norm_A_inv)^(-α)
    a = min(a1, a2)
    b1 = (ϵ*π*(1-α)*(2-α) / (4*sin(α*π)*(3-2α)*norm_A))^(α/(α-1))
    b2 = (2*norm_A)^α
    b = max(b1, b2)
    l = asinh(2*log(a)/(α*π))
    r = asinh(2*log(b)/(α*π))
    return l, r
end


"""
powm_de(A, α, m, l, r)

Computing A^α based on the m-point DE formula on [l,r].
"""
function powm_de(A, α, m, l, r)
    T = typeof(A[1,1])
    F(x) = cosh(x) * exp(α*π/2*sinh(x)) * inv(exp(π/2*sinh(x))*I + A)

    x = LinRange{T}(l, r, m)
    h = (r-l)/(m-1)
    w = fill(h, m)
    w[1], w[m] = h/2, h/2
    n = size(A, 1)
    T1, T2 = zeros(T, n, n), zeros(T, n, n)
    for k = 1:m÷2
        T1 .+= w[k] * F(x[k])
    end
    for k = m:-1:m÷2+1
        T2 .+= w[k] * F(x[k])
    end
    return sin(α*π)/2 * A * (T1 + T2)
end


"""
powm_de(A, α, ϵ; m0=8, use_arpack=false, Ref=nothing)

Computing A^α via adaptive quadrature algorithm based on the DE formula.
The details is in Algorithm 2.
"""
function powm_de(A::Array{TA,2}, α, ϵ; m0=8, use_arpack=false, Exact=nothing, m_max=1000) where {TA}
    A_f64 = convert(Array{Float64,2}, A)
    σ = svdvals(A_f64)
    σ_max, σ_min = maximum(σ), minimum(σ)
    norm_A, norm_A_inv = σ_max, 1/σ_min
    l, r = get_interval(norm_A, norm_A_inv, α, ϵ)
    l, r = TA(l), TA(r)
    π_TA = TA(π)
    F(x) = exp(α*π_TA*sinh(x)/2) * cosh(x) * inv(exp(π_TA*sinh(x)/2)*I + A)

    m = m0
    h = (r-l) / (m-1)
    x = LinRange(l, r, m)
    w = fill(h, m)
    w[1], w[m] = h/2, h/2
    T1, T2 = zero(A), zero(A)
    for k = 1:m÷2
        T1 .+= w[k] * F(x[k])
    end
    for k = m:(-1):m÷2+1
        T2 .+= w[k] * F(x[k])
    end
    T_old = T1 + T2
    T_new = zero(T_old)

    Approx_old = convert(Array{Float64,2}, A*T_old)
    Approx_new = convert(Array{Float64,2}, T_new)
    err_est = (sin(α*π)/2) * opnorm(Approx_new - Approx_old)
    estimates = []
    m_list = [m]
    errors = []
    if !(isnothing(Exact))
        abserr = opnorm(convert(Array{Float64,2}, (sin(α*π)/2)*Approx_old - Exact))
        push!(errors, abserr)
    end

    while err_est > ϵ/2 && m < m_max
        h = h/2
        T1 .= 0
        T2 .= 0
        for k = 1:(m-1)÷2
            T1 .+= h * F(l + (2k-1)*h)
        end
        for k = (m-1):-1:((m-1)÷2+1)
            T2 .+= h * F(l + (2k-1)*h)
        end
        T_new .= T_old/2 + (T1 + T2)

        Approx_new = convert(Array{Float64,2}, A*T_new)
        err_est = (sin(α*π)/2) * opnorm(Approx_new - Approx_old)
        if !(isnothing(Exact))
            abserr = opnorm(convert(Array{Float64,2}, (sin(α*π)/2)*Approx_new - Exact))
            push!(errors, abserr)
        end
        push!(estimates, err_est)
        push!(m_list, 2m-1)

        m = 2m-1
        T_old .= T_new
        Approx_old .= Approx_new
    end
    return errors, estimates, m_list
end
