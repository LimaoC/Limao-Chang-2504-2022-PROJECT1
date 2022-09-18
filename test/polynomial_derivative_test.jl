#############################################################################
#############################################################################
#
# This file contains units tests for polynomial derivative (and product)
# operations
#                                                                               
#############################################################################
#############################################################################

"""
Test derivative (and product) of dense polynomials.
"""
function prod_derivative_test_poly_dense(; N::Int=10^2, seed::Int=0)
    Random.seed!(seed)
    for _ in 1:N
        p1 = rand(PolynomialDense)
        p2 = rand(PolynomialDense)
        p1d = derivative(p1)
        p2d = derivative(p2)
        @assert (p1d * p2) + (p1 * p2d) == derivative(p1 * p2)
    end
    println("prod_derivative_test_poly_dense - PASSED")
end

"""
Test derivative (and product) of sparse polynomials.
"""
function prod_derivative_test_poly_sparse(; N::Int=10^2, seed::Int=0)
    Random.seed!(seed)
    for _ in 1:N
        p1 = rand(PolynomialSparse)
        p2 = rand(PolynomialSparse)
        p1d = derivative(p1)
        p2d = derivative(p2)
        @assert (p1d * p2) + (p1 * p2d) == derivative(p1 * p2)
    end
    println("prod_derivative_test_poly_sparse - PASSED")
end

"""
Test derivative (and product) of sparse big int polynomials.
"""
function prod_derivative_test_poly_sparse_bi(; N::Int=10^2, seed::Int=0)
    Random.seed!(seed)
    for _ in 1:N
        p1 = rand(PolynomialSparseBI)
        p2 = rand(PolynomialSparseBI)
        p1d = derivative(p1)
        p2d = derivative(p2)
        @assert (p1d * p2) + (p1 * p2d) == derivative(p1 * p2)
    end
    println("prod_derivative_test_poly_sparse_bi - PASSED")
end
