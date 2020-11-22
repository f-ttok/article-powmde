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