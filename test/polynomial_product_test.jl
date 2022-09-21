#############################################################################
#############################################################################
#
# This file contains units tests for polynomial product operations
#                                                                               
#############################################################################
#############################################################################

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

    # for _ in 1:N
    #     p_base = PolynomialSparseBI(Term(1, 0))
    #     for _ in 1:N_prods
    #         p = rand(PolynomialSparseBI)
    #         prod = p_base * p
    #         @assert leading(prod) == leading(p_base) * leading(p)
    #         p_base = prod
    #     end
    # end
    println("prod_test_poly_sparse_bi - PASSED")
end
