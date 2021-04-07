using DelimitedFiles
using GSL
using LinearAlgebra
using MatrixMarket
using Printf

cd(@__DIR__)
include("./algorithms/gauss.jl")


"""
    τ(λmax, λmin, α, m)

Compute the scaling parameter τ to compute A^α by using GJ2.
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



function main()
    m_list = 10:10:200
    ϵ = 2.0^(-53)
    Data = Array{Any,2}(undef, length(m_list)+1, 7)
    Data .= ""
    Data[1,1] = "m"
    Data[2:end, 1] .= m_list

    j = 2
    for matname in ["SPD_well", "SPD_ill"]
        println("--- $(matname) ---")
        A = readdlm("matrix/test4_A_$(matname).txt")
        σ = svdvals(A)
        σmax, σmin = maximum(σ), minimum(σ)
        for αx10 in [2, 5, 8]
            print("α = 0.$(αx10): ")
            α = αx10/10
            Exact = readdlm("matrix/test4_Exact$(αx10)_$(matname).txt")
            norm_Exact = norm(Exact)

            Data[1, j] = "$(matname)_$(αx10)_gj2pre"
            for (i, m) in enumerate(m_list)
                print("#")
                c = τ(σmax, σmin, α, m)
                X = powm_gj2(c*A, α, m) / c^α
                Data[i+1, j] = norm(X - Exact) / norm_Exact

            end
            println()
            j += 1
        end
        println()
    end
    writedlm("result/testb1.csv", Data, ',')
end

main()