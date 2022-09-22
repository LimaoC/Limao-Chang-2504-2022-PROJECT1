#############################################################################
#############################################################################
#
# This file contains units tests for testing the sparsity of sparse
# polynomials
#                                                                               
#############################################################################
#############################################################################

function polynomial_sparse_tests()
    poly_sparse_sparsity_construction()
    poly_sparse_bi_sparsity_construction()
    poly_sparse_sparsity_addition()
    poly_sparse_bi_sparsity_addition()
    poly_sparse_sparsity_subtraction()
    poly_sparse_bi_sparsity_subtraction()
end

"""
Test that a PolynomialSparse's terms are correctly constructed
"""
function poly_sparse_sparsity_construction()
    x = x_poly(PolynomialSparse)
    p = x^1000 + 2x + 5
    @assert length(p.terms) == 3
    println("poly_sparse_sparsity_construction - PASSED")
end
function poly_sparse_bi_sparsity_construction()
    x = x_poly(PolynomialSparseBI)
    p = x^1000 + 2x + 5
    @assert length(p.terms) == 3
    println("poly_sparse_bi_sparsity_construction - PASSED")
end

"""
Test that adding zero terms doesn't affect the terms of a PolynomialSparse
"""
function poly_sparse_sparsity_addition(; N::Int=100, seed::Int = 0)
    Random.seed!(seed)
    for _ in 1:N
        p = rand(PolynomialSparse)
        p_terms_before = p.terms
        p += Term(0, degree(p)+1)
        p_terms_after = p.terms
        @assert p_terms_before == p_terms_after
    end
    println("poly_sparse_sparsity_addition - PASSED")
end
function poly_sparse_bi_sparsity_addition(; N::Int=100, seed::Int = 0)
    Random.seed!(seed)
    for _ in 1:N
        p = rand(PolynomialSparseBI)
        p_terms_before = p.terms
        p += Term(0, degree(p)+1)
        p_terms_after = p.terms
        @assert p_terms_before == p_terms_after
    end
    println("poly_sparse_bi_sparsity_addition - PASSED")
end

"""
Test that zeroing out terms (via subtraction) in a PolynomialSparse removes the term
"""
function poly_sparse_sparsity_subtraction(; N::Int=100, seed::Int = 0)
    Random.seed!(seed)
    for _ in 1:N
        p = rand(PolynomialSparse)
        iszero(p) && break  # zero polynomial will break this test
        p_num_terms_before = length(p.terms)
        # zero out leading term
        p -= leading(p)
        p_num_terms_after = length(p.terms)
        @assert p_num_terms_before > p_num_terms_after
    end
    println("poly_sparse_sparsity_addition - PASSED")
end
function poly_sparse_bi_sparsity_subtraction(; N::Int=100, seed::Int = 0)
    Random.seed!(seed)
    for _ in 1:N
        p = rand(PolynomialSparseBI)
        iszero(p) && break  # zero polynomial will break this test
        p_num_terms_before = length(p.terms)
        # zero out leading term
        p -= leading(p)
        p_num_terms_after = length(p.terms)
        @assert p_num_terms_before > p_num_terms_after
    end
    println("poly_sparse_bi_sparsity_addition - PASSED")
end
