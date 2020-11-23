using LinearAlgebra
using Printf


"""
    sqrtm_newtom(A; tol=eps(), maxiter=50)

Computing the (principal) square root of A
based on the scaled Denman--Beavers iteration [Eq(6.29), 1].

[1]: Functions of matrices: Theory and computation. SIAM, Philadelphia, PA, 2008.
"""
function sqrtm_newton(A; tol=eps(), maxiter=50)
    println("Computing the square root of A.")
    n = size(A, 1)
    X = copy(A)
    M = copy(A)
    M_inv = zero(A)
    norm_A = norm(A)

    iter = 0
    relres = norm(X^2 - A) / norm_A
    while iter < maxiter && relres > tol
        iter += 1
        μ = relres > 1e-2 ? (det(M))^(-1/(2n)) : 1.0
        M_inv .= inv(M)
        X .= μ * X * (I + μ^(-2)*M_inv) / 2
        M .= (I + (μ^2*M + μ^(-2)*M_inv)/2) / 2
        relres = norm(X^2 - A) / norm_A
    end
    @printf("\titer: %d, relres: %.3e\n", iter, relres)
    return X
end


"""
prootm_newton(A; tol=eps(), maxiter=50)

Computing the (principal) pth root of A
based on the inverse Newton iteration [Alg. 7.14, 1]

[1]: Functions of matrices: Theory and computation. SIAM, Philadelphia, PA, 2008.
"""
function prootm_newton(A::Array{TA,2}; tol=eps(), maxiter=50) where {TA}
    B = sqrtm_newton(A, tol=tol, maxiter=maxiter)
    B .= sqrtm_newton(B, tol=tol, maxiter=maxiter)
    θ = opnorm(B, 1)
    c = (θ/sqrt(2))^(1/p)

    X = diagm(ones(TA, n)) / c
    M = B / c^p
    T = zero(B)

    println("Computing the $(-p)th root of A^(1/4)")
    iter = 0
    relres = norm(B*X^p - I) / sqrt(n)
    while iter < maxiter && relres > tol
        iter += 1
        T .= ((p+1)*I - M) / p
        X .= X * T
        M .= T^p * M
        relres = norm(B*X^p - I) / sqrt(n)
    end
    @printf("\titer: %d, relres: %.3e\n", iter, relres)

    X .= (inv(X))^4
    @printf("relres (total): %.3e\n", norm(A - X^p) / norm(A))
    return X
end