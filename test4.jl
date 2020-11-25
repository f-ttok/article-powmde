using LinearAlgebra
using DelimitedFiles
using Printf

cd(@__DIR__)

function generate_exact_solutions()
    for matname in 
end

function compute_refsol()
    println("compute reference solutions")
    matnamelist = ["SPD_well", "SPD_ill", "NS_well", "NS_ill"]
    for matname in matnamelist
        println("$(matname)")
        A_f64 = readdlm("matrices/test3_A_$(matname).txt")
        A = convert(Array{BigFloat, 2}, A_f64)
        n = size(A,1)
        norm_A = norm(A)
        B = (A / norm_A)

        maxiter = 100
        X = diagm(0 => ones(BigFloat, n))
        N = copy(B)
        relres = norm(X^2 - B)
        maxiter = 100
        iter = 0
        tol = 1e-50

        println("Square root")
        while iter < maxiter && relres > tol
            @printf("%4d | %.3e\n", iter, relres)
            iter += 1
            T = (I+N)/2
            X = X*T
            N = (T^2) \ N
            relres = norm(X^2 - B)
        end
        println()

        C = copy(X)
        X = diagm(0 => ones(BigFloat, n))
        N = copy(C)
        norm_C = norm(C)
        relres = norm(X^5 - C) / norm_C
        maxiter = 100
        iter = 0
        tol = 1e-50
        println("Fifth root")
        while iter < maxiter && relres > tol
            @printf("%4d | %.3e\n", iter, relres)
            iter += 1
            T = (4I+N)/5
            X = X*T
            N = (T^5) \ N
            relres = norm(X^5 - C) / norm_C
        end

        Exact = X * norm_A^(1/10)
        Y = copy(Exact)
        for i = 2:8
            Y = Y * Exact
            if i in [2,5,8]
                writedlm("matrices/test3_Exact_$(matname)_$(i).txt", convert(Array{Float64,2}, Y))
            end
        end
    end
end


function main()
    for αx10 in [2, 5, 8]
        matnamelist = ["SPD_well", "SPD_ill", "NS_well", "NS_ill"]
        α = αx10/10
        for matname in matnamelist
            println("--- $(matname) ---")
            A = readdlm("matrices/test3_A_$(matname).txt")
            Ref = readdlm("matrices/test3_Exact_$(matname)_$(αx10).txt")
            norm_Ref = norm(Ref)
            ρ = maximum(abs.(eigvals(A)))

            l, r = get_lr(A, α, 2.0^(-53)*ρ^α, use_arpack=false)

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
            writedlm("results/test3_$(matname)_$(αx10).txt", Errlist)
            println()
        end
    end
end