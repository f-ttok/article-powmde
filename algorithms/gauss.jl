using LinearAlgebra
using FastGaussQuadrature

function powm_gj1(A, α, m)
    u, w = gaussjacobi(m, 1/α-2, 0)
    F(u) = inv((1+u)^(1/α)*I + (1-u)^(1/α)*A)

    G = zero(A)
    for k = 1:m
        G .= G .+ w[k] * F(u[k])
    end
    G .= 2*sin(α*π) / (α*π) .* (A * G)
    return G
end


function powm_gj2(A, α, m)
    v, w = gaussjacobi(m, α-1, -α)
    F(v) = inv((1-v)*I + (1+v)*A)

    G = zero(A)
    for k = 1:m
        G .= G .+ w[k] * F(v[k])
    end
    G .= 2*sin(α*π) / π .* (A * G)
    return G
end
