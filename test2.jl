using Arpack
using DelimitedFiles
using LinearAlgebra
using MatrixMarket
using Printf
using SparseArrays


include("./algorithms/newton.jl")
include("./algorithms/de.jl")


function generate_reference_solutions()
    for matname in ["ex5", "pores_1"]
        A_f64 = collect(mmread("matrix/$(matname).mtx"))
        if matname == "pores_1"
            A_f64 = -1.0 * A_f64
        end
        writedlm("matrix/test2_A_$(matname).txt", A_f64)
        A = convert(Array{BigFloat,2}, A_f64)
        n = size(A,1)
        tol = sqrt(n) * eps(BigFloat)/2
        # Exact2 = prootm_newton(A)
        # Exact5 = sqrtm_newton(A, tol=tol, maxiter=100)
        # Exact5_f64 = convert(Array{Float64,2}, Exact)
        # writedlm("matrix/test1_Exact_$(matname).txt", Exact)
    end
end
