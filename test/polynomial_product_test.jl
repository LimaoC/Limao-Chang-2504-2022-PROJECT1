#############################################################################
#############################################################################
#
# This file contains units tests for polynomial product operations
#                                                                               
#############################################################################
#############################################################################

"""
Execute all polynomial product tests in this file.
"""
function polynomial_product_tests()
    @time prod_test_poly_dense()
    @time prod_test_poly_sparse()
    @time prod_test_poly_sparse_bi()
    @time prod_test_poly_mod_p()
end

"""
Test product of dense polynomials.
"""
function prod_test_poly_dense(; N::Int=10^3, N_prods::Int=20, seed::Int=0)
    Random.seed!(seed)
    for _ in 1:N
        p1 = rand(PolynomialDense)
        p2 = rand(PolynomialDense)
        prod = p1 * p2
        @assert leading(prod) == leading(p1) * leading(p2)
    end

    for _ in 1:N
        p_base = PolynomialDense(Term(1, 0))
        for _ in 1:N_prods
            p = rand(PolynomialDense)
            prod = p_base * p
            @assert leading(prod) == leading(p_base) * leading(p)
            p_base = prod
        end
    end
    println("prod_test_poly_dense - PASSED")
end

"""
Test product of sparse polynomials.
"""
function prod_test_poly_sparse(; N::Int=10^3, N_prods::Int=20, seed::Int=0)
    Random.seed!(seed)
    for _ in 1:N
        p1 = rand(PolynomialSparse)
        p2 = rand(PolynomialSparse)
        prod = p1 * p2
        @assert leading(prod) == leading(p1) * leading(p2)
    end

    for _ in 1:N
        p_base = PolynomialSparse(Term(1, 0))
        for _ in 1:N_prods
            p = rand(PolynomialSparse)
            prod = p_base * p
            @assert leading(prod) == leading(p_base) * leading(p)
            p_base = prod
        end
    end
    println("prod_test_poly_sparse - PASSED")
end

"""
Test product of sparse big int polynomials.
"""
function prod_test_poly_sparse_bi(; N::Int=100, N_prods::Int=20, seed::Int=0)
    Random.seed!(seed)
    for _ in 1:N
        p1 = rand(PolynomialSparseBI)
        p2 = rand(PolynomialSparseBI)
        prod = p1 * p2
        @assert leading(prod) == leading(p1) * leading(p2)
    end

    for _ in 1:N
        p_base = PolynomialSparseBI(Term(1, 0))
        for _ in 1:N_prods
            p = rand(PolynomialSparseBI)
            prod = p_base * p
            @assert leading(prod) == leading(p_base) * leading(p)
            p_base = prod
        end
    end
    println("prod_test_poly_sparse_bi - PASSED")
end

"""
Test product of polynomials modulo some prime.
"""
function prod_test_poly_mod_p(; N::Int=100, N_prods::Int=20, seed::Int=0)
    Random.seed!(seed)
    for _ in 1:N
        # make sure we are using the same prime for p1 and p2
        rand_prime = prime(rand(1:100))
        p1 = rand(PolynomialModP, p=rand_prime)
        p2 = rand(PolynomialModP, p=rand_prime)
        prod = p1 * p2
        @assert leading(prod) == mod(leading(p1) * leading(p2), rand_prime)
    end

    for _ in 1:N
        rand_prime = prime(rand(1:100))
        p_base = PolynomialModP(PolynomialSparse(Term(1, 0)), rand_prime)
        for _ in 1:N_prods
            p = rand(PolynomialModP, p=rand_prime)
            prod = p_base * p
            @assert (leading(prod.polynomial) ==
                mod(leading(p_base.polynomial) * leading(p.polynomial), rand_prime))
            p_base = prod
        end
    end
    println("prod_test_poly_mod_p - PASSED")
end
