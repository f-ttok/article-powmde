using Arpack
using DelimitedFiles
using LinearAlgebra
using MatrixMarket
using Printf
using SparseArrays

cd(@__DIR__)

include("./algorithms/newton.jl")
include("./algorithms/de.jl")


function generate_exact_solutions()
    for matname in ["ex5", "pores_1"]
        println("--- $(matname) ---")
        A_f64 = collect(mmread("matrix/$(matname).mtx"))
        if matname == "pores_1"
            A_f64 = -1.0 * A_f64
        end
        writedlm("matrix/test2_A_$(matname).txt", A_f64)
        A = convert(Array{BigFloat,2}, A_f64)
        n = size(A,1)
        tol = sqrt(n) * eps(BigFloat)/2
        println("Computing A^0.2")
        Exact2 = prootm_newton(A, 5, tol=tol)
        println("Computing A^0.5")
        Exact5 = sqrtm_newton(A, tol=tol)
        Exact8 = Exact2^4
        writedlm("matrix/test2_Exact2_$(matname).txt", Exact2)
        writedlm("matrix/test2_Exact5_$(matname).txt", Exact5)
        writedlm("matrix/test2_Exact8_$(matname).txt", Exact8)
        println()
    end
end


function main()
    Data = Array{Any, 2}(undef, 9, 25)
    Data .= ""
    Data[1,1] = "m"
    m_list = [7*(2^k)+1 for k=0:7]
    Data[2:end,1] .= m_list

    j = 2
    for matname in ["ex5", "pores_1"]
        println("--- $(matname) ---")
        A = readdlm("matrix/test2_A_$(matname).txt")
        λ = eigvals(A)
        ρ = maximum(abs.(λ))
        σ = svdvals(A)
        σmax, σmin = BigFloat(maximum(σ)), BigFloat(minimum(σ))
        c = 1 / sqrt(σmax*σmin)
        A_tilde = convert(Array{BigFloat,2}, c*A)

        for αx10 in [2, 5, 8]
            println("\tα = 0.$(αx10)")
            α = αx10/10
            Exact = readdlm("matrix/test2_Exact$(αx10)_$(matname).txt", '\t', BigFloat)
            Exact_tilde = c^α * Exact
            norm_Exact = opnorm(convert(Array{Float64,2}, Exact))
            for log10ϵ in [-7,-14]
                ϵ = 10.0^log10ϵ
                ϵ_tilde = ϵ * c^α * ρ^α
                Data[1,j] = "$(matname)_$(αx10)_$(log10ϵ)_est"
                Data[1,j+1] = "$(matname)_$(αx10)_$(log10ϵ)_err"
                errors, estimates, m_list = powm_de(A_tilde, α, ϵ_tilde, m0=8, Exact=Exact_tilde)
                n = length(estimates)
                Data[3:n+2, j] .= convert(Array{Float64,1}, estimates / c^α / norm_Exact)
                n = length(errors)
                Data[2:n+1, j+1] .= convert(Array{Float64,1}, errors / c^α / norm_Exact)
                j += 2
            end
        end
    end
    writedlm("result/test2.csv", Data, ',')
end
