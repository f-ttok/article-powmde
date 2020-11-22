using LinearAlgebra
using Printf


"""
    sqrtm_newtom(A; tol=eps(), maxiter=50)

Computing the (principal) square root of A
based on the scaled Denman--Beavers iteration [6.29, 1].

[1]: Functions of matrices: Theory and computation. SIAM, Philadelphia, PA, 2008.
"""
function sqrtm_newton(A; tol=eps(), maxiter=50)
    println("Computing the square root of A.")
    n = size(A, 1)
    X = copy(A)
    M = copy(A)
    norm_A = norm(A)

    iter = 0
    relres = norm(X^2 - A) / norm_A
    while iter < maxiter && relres > tol
        iter += 1
        μ = relres > 1e-2 ? (det(M))^(-1/(2n)) : 1.0
        M_inv = inv(M)
        X = μ * X * (I + μ^(-2)*M_inv) / 2
        M = (I + (μ^2*M + μ^(-2)*M_inv)/2) / 2
        relres = norm(X^2 - A) / norm_A
    end
    @printf("iter: %d, relres: %.3e\n", iter, relres)
    return X
end