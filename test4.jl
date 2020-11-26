using LinearAlgebra
using DelimitedFiles
using MatrixMarket
using Printf

cd(@__DIR__)
include("./algorithms/de.jl")
include("./algorithms/gauss.jl")
include("./algorithms/newton.jl")

function generate_exact_solutions()
    for matname in ["SPD_well", "SPD_ill", "NS_well", "NS_ill"]
        println("--- $(matname) ---")
        A_f64 = readdlm("matrix/test4_A_$(matname).txt")
        A = convert(Array{BigFloat,2}, A_f64)
        n = size(A,1)
        tol = sqrt(n) * eps(BigFloat)/2
        println("Computing A^0.2")
        Exact2 = prootm_newton(A, 5, tol=tol, maxiter=50)
        println("Computing A^0.5")
        Exact5 = sqrtm_newton(A, tol=tol, maxiter=50)
        Exact8 = Exact2^4
        writedlm("matrix/test4_Exact2_$(matname).txt", convert(Array{Float64,2}, Exact2))
        writedlm("matrix/test4_Exact5_$(matname).txt", convert(Array{Float64,2}, Exact5))
        writedlm("matrix/test4_Exact8_$(matname).txt", convert(Array{Float64,2}, Exact8))
        println()
    end
end


function main()
    m_list = 10:10:200
    ϵ = 2.0^(-53)
    Data = Array{Any,2}(undef, length(m_list)+1, 37)
    Data .= ""
    Data[1,1] = "m"
    Data[2:end, 1] .= m_list

    j = 2
    for matname in ["SPD_well", "SPD_ill", "NS_well", "NS_ill"]
        println("--- $(matname) ---")
        A = readdlm("matrix/test4_A_$(matname).txt")
        σ = svdvals(A)
        σmax, σmin = maximum(σ), minimum(σ)
        c = 1 / sqrt(σmax*σmin)
        A_tilde = c*A
        norm_A_tilde = c * σmax
        norm_A_tilde_inv = c / σmin
        ρ = maximum(abs.(eigvals(A)))
        for αx10 in [2, 5, 8]
            print("α = 0.$(αx10): ")
            α = αx10/10
            Exact = readdlm("matrix/test4_Exact$(αx10)_$(matname).txt")
            norm_Exact = norm(Exact)

            ϵ_tilde = ϵ * c^α * ρ^α
            l, r = get_interval(norm_A_tilde, norm_A_tilde_inv, α, ϵ_tilde)
            Data[1, j] = "$(matname)_$(αx10)_gj1"
            Data[1, j+1] = "$(matname)_$(αx10)_gj2"
            Data[1, j+2] = "$(matname)_$(αx10)_de"
            for (i, m) in enumerate(m_list)
                print("#")
                X = powm_gj1(A_tilde, α, m) / c^α
                Data[i+1, j] = norm(X - Exact) / norm_Exact

                X = powm_gj2(A_tilde, α, m) / c^α
                Data[i+1, j+1] = norm(X - Exact) / norm_Exact

                X = powm_de(A_tilde, α, m, l, r) / c^α
                Data[i+1, j+2] = norm(X - Exact) / norm_Exact
            end
            println()
            j += 3
        end
        println()
    end
    writedlm("result/test4.csv", Data, ',')
end
