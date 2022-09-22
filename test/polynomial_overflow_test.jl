#############################################################################
#############################################################################
#
# This file contains units tests for polynomial operations that (could) overflow
#                                                                               
#############################################################################
#############################################################################

"""
Execute all polynomial overflow tests in this file.
"""
function polynomial_overflow_tests()
    # the PolynomialSparse tests are meant to fail, so they are commented out to avoid the
    # test erroring and exiting

    # poly_sparse_overflow()
    @time poly_sparse_bi_overflow()
    # poly_sparse_underflow()
    @time poly_sparse_bi_underflow()
    # poly_sparse_prod_overflow()
    @time poly_sparse_bi_prod_overflow()
end

"""
Test whether a polynomial overflows. (note: this test should fail for sparse polynomials)
"""
function poly_sparse_overflow(; N::Int = 128)
    p = x_poly(PolynomialSparse)
    for _ in 1:N
        @assert leading(p * 2) > leading(p)
        p *= 2
    end
    println("poly_sparse_overflow - PASSED")
end

"""
Test whether a polynomial overflows.
"""
function poly_sparse_bi_overflow(; N::Int = 128)
    p = x_poly(PolynomialSparseBI)
    for _ in 1:N
        @assert leading(p * 2) > leading(p)
        p *= 2
    end
    println("poly_sparse_bi_overflow - PASSED")
end

"""
Test whether a polynomial underflows. (note: this test should fail for sparse polynomials)
"""
function poly_sparse_underflow(; N::Int = 128)
    p = -x_poly(PolynomialSparse)
    for _ in 1:N
        @assert leading(p * 2) < leading(p)
        p *= 2
    end
    println("poly_sparse_underflow - PASSED")
end

"""
Test whether a polynomial overflows.
"""
function poly_sparse_bi_underflow(; N::Int = 128)
    p = -x_poly(PolynomialSparseBI)
    for _ in 1:N
        @assert leading(p * 2) < leading(p)
        p *= 2
    end
    println("poly_sparse_bi_underflow - PASSED")
end

"""
Test whether the leading term of a product of two polynomials has degree equal to the sum
of the leading terms of the two polynomials. (note: this test should fail for sparse
polynomials)
"""
function poly_sparse_prod_overflow()
    x = x_poly(PolynomialSparse)
    p1 = (2^32)x^2 + 5x + 3
    p2 = (2^40)x^3 + 9x^2 + 2x
    @assert leading(p1 * p2).degree == leading(p1).degree + leading(p2).degree
    println("poly_sparse_prod_overflow - PASSED")
end

"""
Test whether the leading term of a product of two polynomials has degree equal to the sum
of the leading terms of the two polynomials.
"""
function poly_sparse_bi_prod_overflow()
    x = x_poly(PolynomialSparseBI)
    p1 = (2^32)x^2 + 5x + 3
    p2 = (2^40)x^3 + 9x^2 + 2x
    @assert leading(p1 * p2).degree == leading(p1).degree + leading(p2).degree
    println("poly_sparse_bi_prod_overflow - PASSED")
end