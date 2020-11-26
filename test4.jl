using LinearAlgebra
using DelimitedFiles
using Printf

cd(@__DIR__)
include("./algorithms/de.jl")
include("./algorithms/newton.jl")

function generate_exact_solutions()
    for matname in ["SPD_well", "SPD_ill", "NS_well", "NS_ill"]
        println("--- $(matname) ---")
        A_f64 = collect(mmread("matrix/$(matname).mtx"))
        writedlm("matrix/test4_A_$(matname).txt", A_f64)
        A = convert(Array{BigFloat,2}, A_f64)
        n = size(A,1)
        tol = sqrt(n) * eps(BigFloat)/2
        println("Computing A^0.2")
        Exact2 = prootm_newton(A, 5, tol=tol)
        println("Computing A^0.5")
        Exact5 = sqrtm_newton(A, tol=tol)
        Exact8 = Exact2^4
        writedlm("matrix/test4_Exact2_$(matname).txt", convert(Array{Float64,2}, Exact2))
        writedlm("matrix/test4_Exact5_$(matname).txt", convert(Array{Float64,2}, Exact5))
        writedlm("matrix/test4_Exact8_$(matname).txt", convert(Array{Float64,2}, Exact8))
        println()
    end
end


function main()
    m_list = 10:10:200
    Data = Array{Any,2}(undef, length(m_list)+1, 37)
    Data
    for αx10 in [2, 5, 8]
        matnamelist = ["SPD_well", "SPD_ill", "NS_well", "NS_ill"]
        α = αx10/10
        for matname in matnamelist
            println("--- $(matname) ---")
            A = readdlm("matrix/test4_A_$(matname).txt")
            Ref = readdlm("matrix/test4_Exact$(matname)_$(αx10).txt")
            norm_Ref = norm(Ref)
            ρ = maximum(abs.(eigvals(A)))

            l, r = get_interval

            mlist = 10:10:200
            Errlist = zeros(length(mlist), 4)
            Errlist[:,1] .= mlist
            for (i, m) in enumerate(mlist)
                print("$(m) ")
                X = powm_gj1(A, α, m)
                Errlist[i,2] = norm(X - Ref) / norm_Ref

                X = powm_gj2(A, α, m)
                Errlist[i,3] = norm(X - Ref) / norm_Ref

                X = powm_de(A, α, m, l, r)
                Errlist[i,4] = norm(X - Ref) / norm_Ref
            end
            writedlm("result/test4_$(matname)_$(αx10).txt", Errlist)
            println()
        end
    end
end