using DelimitedFiles
using FastGaussQuadrature
using LinearAlgebra
using MatrixMarket
using Arpack
using SparseArrays
using Printf

cd(@__DIR__)

include("./algorithms/de.jl")
include("./algorithms/newton.jl")

function generate_exact_solutions()
    for matname in ["lund_b", "bcsstk04", "nos4"]
        A_f64 = collect(mmread("matrix/$(matname).mtx"))
        writedlm("matrix/test3_A_$(matname).txt", A_f64)
        A = convert(Array{BigFloat,2}, A_f64)
        n = size(A,1)
        tol = sqrt(n) * eps(BigFloat)/2
        Exact = sqrtm_newton(A, tol=1e-20, maxiter=50)
        Exact = convert(Array{Float64,2}, Exact)
        writedlm("matrix/test3_Exact_$(matname).txt", Exact)
    end
end


function main()
    m_list = 10:10:200
    Data = Array{Any,2}(undef, length(m_list)+1, 7)
    Data .= ""
    Data[1,1] = "m"
    Data[2:end,1] .= m_list
    α = 0.5
    j = 2
    ϵ = 2.0^(-53)
    for matname in ["lund_b", "bcsstk04", "nos4"]
        println("----- $(matname) -----")
        A = readdlm("matrix/test3_A_$(matname).txt")
        λ = eigvals(A)
        λmax, λmin = maximum(λ), minimum(λ)
        c = 1 / sqrt(λmax*λmin)
        A_tilde = c*A
        λmax_tilde = c*λmax
        A1_tilde = ones(1,1) * λmax_tilde
        Exact_tilde = c^α * readdlm("matrix/test3_Exact_$(matname).txt")
        Exact1_tilde = ones(1,1) * sqrt(λmax_tilde)

        μ = log(λmax_tilde)^2 + π^2 + 1
        d0 = asin(sqrt((μ - sqrt(μ^2 - 4π^2)) / 2))
        ϵ_tilde = ϵ * c^α * λmax_tilde^α
        l, r = get_interval(λmax_tilde, λmax_tilde, α, ϵ_tilde)
        @show 2π*d0/(r-l)

        Data[1,j] = "$(matname)_A"
        Data[1,j+1] = "$(matname)_λ"
        for (i, m) in enumerate(m_list)
            print("$(m) ")
            X = powm_de(A_tilde, α, m, l, r)
            Data[i+1,j] = opnorm(X - Exact_tilde) / λmax_tilde^α
            X = powm_de(A1_tilde, α, m, l, r)
            Data[i+1,j+1] = opnorm(X - Exact1_tilde) / λmax_tilde^α
        end
        j += 2
        println()
    end
    writedlm("result/test3.csv", Data, ',')
end