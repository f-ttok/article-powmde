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
        A_f64 = collect(mmread("matrix/$(matname).mtx"))
        if matname == "pores_1"
            A_f64 = -1.0 * A_f64
        end
        A = convert(Array{BigFloat,2}, A_f64)
        n = size(A,1)
        tol = sqrt(n) * eps(BigFloat)/2
        Exact = sqrtm_newton(A, tol=tol, maxiter=100)
        Exact = convert(Array{Float64,2}, Exact)
        writedlm("matrix/test1_A_$(matname).txt", A_f64)
        writedlm("matrix/test1_Exact_$(matname).txt", Exact)
    end
end


function main()
    m_list = 10:10:200
    Data = Array{Any,2}(undef, length(m_list)+1, 15)
    Data .= ""
    Data[1,1] = "m"
    Data[2:end,1] .= m_list

    Data_interval = Array{Any,2}(undef, 5, 6)
    Data_interval .= ""
    Data_interval[1,:] = ["matrix", "log10ϵ", "l", "r_svdvals", "r_arpack", "r_arpackmod"]

    α = 0.5
    Δ = 1e-3  # (relative) tolerance for eigenvalues/singular values
    ncv = 3  # dimensions of the Krylov subspace for computing
             # eigenvalues/singular values

    i_interval = 2
    j = 2
    for matname in ["ex5", "pores_1"]
        println("----- $(matname) -----")
        A_f64 = readdlm("matrix/test1_A_$(matname).txt")
        σ = svdvals(A_f64)
        σmax, σmin = maximum(σ), minimum(σ)
        c = 1 / sqrt(σmax * σmin)
        ρ = maximum(abs.(eigvals(A_f64)))
        A = convert(Array{BigFloat,2}, A_f64)
        A_f64_tilde = c * A_f64
        A_tilde = c * A
        Exact = readdlm("matrix/test1_Exact_$(matname).txt")

        X, X_f64 = zero(A), zero(A_f64)

        for log10ϵ in [-7, -14]
            ϵ = 10.0 ^ log10ϵ
            @printf("tol: 10^%d\n", log10ϵ)
            ϵ_tilde = ϵ * ρ^α * c^α

            println("\tcomputing intervals...")
            # Computing [l,r] based on svdvals
            norm_A, norm_A_inv = σmax*c, σmax*c
            l, r = get_interval(norm_A, norm_A_inv, α, ϵ_tilde)

            # Computing [l,r] based on Arpack.jl
            if issymmetric(A)
                B = sparse(A_f64_tilde)
                norm_A = eigs(B, nev=1, tol=Δ, which=:LM, ncv=ncv)[1][1]
                norm_A_inv = 1 / eigs(B, nev=1, tol=Δ, which=:SM, ncv=ncv)[1][1]
            else
                B = sparse(A_f64_tilde)
                BtB = B'*B
                norm_A = sqrt(real(eigs(BtB, nev=1, tol=Δ, which=:LM, ncv=ncv)[1][1]))
                norm_A_inv = 1 / sqrt(real(eigs(BtB, nev=1, tol=Δ, which=:SM, ncv=ncv)[1][1]))
            end
            l_ap, r_ap = get_interval(norm_A, norm_A_inv, α, ϵ_tilde)

            # Computing [l,r] based on Arpack.jl with the modified tolerance
            l_apm, r_apm = get_interval(norm_A, norm_A_inv, α, ϵ_tilde/(1 + 1/(1-Δ)))

            Data_interval[i_interval,:] = [matname, log10ϵ, l, r, r_ap, r_apm]
            i_interval += 1

            Data[1,j] = "$(matname)_$(log10ϵ)_svdvals"
            Data[1,j+1] = "$(matname)_$(log10ϵ)_Arpack"
            Data[1,j+2] = "$(matname)_$(log10ϵ)_Arpackmod"

            print("\tcomputing A^0.5 ")
            for (i, m) in enumerate(m_list)
                print("#")
                X .= powm_de(A_tilde, α, m, l, r) / c^α
                X_f64 .= convert(Array{Float64,2}, X)
                Data[i+1, j] = opnorm(X_f64 - Exact) / norm(Exact)

                X .= powm_de(A_tilde, α, m, l_ap, r_ap)
                X_f64 = convert(Array{Float64,2}, X) / c^α
                Data[i+1, j+1] = opnorm(X_f64 - Exact) / norm(Exact)

                X .= powm_de(A_tilde, α, m, l_apm, r_apm)
                X_f64 .= convert(Array{Float64,2}, X) / c^α
                Data[i+1, j+2] = opnorm(X_f64 - Exact) / norm(Exact)
            end
            j += 3
            println()
        end

        Data[1,j] = "$(matname)_wide"
        print("computing A^0.5 with the wide interval\n\t")
        for (i,m) in enumerate(m_list)
            print("#")
            X .= powm_de(A_tilde, α, m, -6.0, 6.0)
            X_f64 .= convert(Array{Float64,2}, X) / c^α
            Data[i+1, j] = opnorm(X_f64 - Exact) / norm(Exact)
        end
        j += 1
        println()
    end
    writedlm("result/test1_interval.csv", Data_interval, ',')
    writedlm("result/test1.csv", Data, ',')
end