#############################################################################
#############################################################################
#
# This file contains units tests for extended euclid algorithm polynomial
# operations
#                                                                               
#############################################################################
#############################################################################

"""
Executes all extended elucid algorithm polynomial tests in this file.
"""
function polynomial_ext_euclid_tests()
    ext_euclid_test_poly_dense()
    ext_euclid_test_poly_sparse()
    ext_euclid_test_poly_sparse_bi()
    ext_euclid_test_poly_mod_p()
end

"""
Test the extended euclid algorithm for dense polynomials modulo p.
"""
function ext_euclid_test_poly_dense(; prime::Int=101, N::Int=10^3, seed::Int=0)
    Random.seed!(seed)
    for _ in 1:N
        p1 = rand(PolynomialDense)
        p2 = rand(PolynomialDense)
        g, s, t = extended_euclid_alg(p1, p2, prime)
        @assert mod(s * p1 + t * p2 - g, prime) == 0
    end
    println("ext_euclid_test_poly_dense - PASSED")
end

"""
Test the extended euclid algorithm for sparse polynomials modulo p.
"""
function ext_euclid_test_poly_sparse(; prime::Int=101, N::Int=10^3, seed::Int=0)
    Random.seed!(seed)
    for _ in 1:N
        p1 = rand(PolynomialSparse)
        p2 = rand(PolynomialSparse)
        g, s, t = extended_euclid_alg(p1, p2, prime)
        @assert mod(s * p1 + t * p2 - g, prime) == 0
    end
    println("ext_euclid_test_poly_sparse - PASSED")
end

"""
Test the extended euclid algorithm for sparse big int polynomials modulo p.
"""
function ext_euclid_test_poly_sparse_bi(; prime::Int=101, N::Int=10^3, seed::Int=0)
    Random.seed!(seed)
    for _ in 1:N
        p1 = rand(PolynomialSparse)
        p2 = rand(PolynomialSparse)
        g, s, t = extended_euclid_alg(p1, p2, prime)
        @assert mod(s * p1 + t * p2 - g, prime) == 0
    end
    println("ext_euclid_test_poly_sparse_bi - PASSED")
end

"""
Test the extended euclid algorithm for polynomials modulo p.
"""
function ext_euclid_test_poly_mod_p(; prime::Int=101, N::Int=10^3, seed::Int=0)
    Random.seed!(seed)
    for _ in 1:N
        p1 = rand(PolynomialModP, p=prime)
        p2 = rand(PolynomialModP, p=prime)
        g, s, t = extended_euclid_alg(p1, p2, prime)
        g = PolynomialModP(g, prime)
        s = PolynomialModP(s, prime)
        t = PolynomialModP(t, prime)
        # @assert mod(s * p1 + t * p2 - g, prime) == 0
        @assert s * p1 + t * p2 - g == 0
    end
    println("ext_euclid_test_poly_mod_p - PASSED")
end