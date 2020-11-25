using DelimitedFiles
using LinearAlgebra
using MatrixMarket
using SparseArrays

"""
    matrices_test4()

Generate test matrices for test4.
In our tests, the parameter c is set as follows:

----- NS_well ------

(c, cond(A)) = (0.08633279800415047, 100.00002382542779)
------ NS_ill ------

(c, cond(A)) = (0.3028873801231384, 1.0000005120888444e7)
"""
function matrices_test4()
    println("generate matrices for test4")
    n = 100
    # R = rand(n,n)
    # writedlm("matrix/test4_R.txt", R)
    R = readdlm("matrix/test4_R.txt")

    Q = qr(R).Q
    for (κ, matname) in zip([2, 7], ["SPD_well", "SPD_ill"])
        d = 10.0 .^ LinRange(-κ/2, κ/2, n)
        D = diagm(0 => d)
        A = Q * D * Q'
        writedlm("matrix/test4_A_$(matname).txt", A)
    end

    for (κ, matname) in zip([1e2, 1e7], ["NS_well", "NS_ill"])
        println("------ $(matname) ------")
        is_converged = false
        c = 1.0
        a = 1.0
        b = 1e-16
        iter = 0

        while !(is_converged)
            iter += 1
            c = (a+b)/2
            A = exp(c*R)
            k = cond(A)
            # @printf("%5d, c = %.3e, κ(A) = %.3e\n", iter, c, k)
            if abs(k - κ)/κ < 1e-6
                is_converged = true
            elseif k - κ > 0
                a = c
            else
                b = c
            end

            if iter > 100
                is_converged = true
            end
        end
        println()
        @show c, cond(A)
        σ = svdvals(A)
        A = A / sqrt(maximum(σ)*minimum(σ))

        writedlm("matrix/test4_A_$(matname).txt", A)
    end

end


function matrices_test5()
    println("generate a matrix for test5")
    n = 200
    d0 = fill(2, n)
    d1 = fill(-1, n-1)
    L = spdiagm(-1=>d1, 0=>d0, 1=>d1)
    In = Diagonal(ones(n))
    A = kron(L, In) + kron(In, L)

    mmwrite("matrix/poisson200.mtx", A)
end


matrices_test4()
matrices_test5()