#############################################################################
#############################################################################
#
# This file contains units tests for PolynomialModP polynomials
#                                                                               
#############################################################################
#############################################################################

"""
Executes all the unit tests for PolynomialModPs in this file.
"""
function polynomial_mod_p_tests()
    poly_mod_p_invalid_sum()
    poly_mod_p_invalid_subtract()
    poly_mod_p_invalid_prod()
    poly_mod_p_invalid_p()
    poly_mod_p_same_poly_diff_prime()
end

"""
Test adding two PolynomialModPs that have different prime p's.
"""
function poly_mod_p_invalid_sum(; seed::Int = 0)
    Random.seed!(seed)
    # try and add two PolynomialModPs that have different p's
    p1 = PolynomialModP(rand(PolynomialSparse), 17)
    p2 = PolynomialModP(rand(PolynomialSparse), 5)
    try
        p1 + p2
    catch AssertionError
        println("poly_mod_p_invalid_sum - PASSED")
        return
    end
    throw("PolynomialModP's that have different p's shouldn't be able to be added " *
          "together")
end

"""
Test subtracting a PolynomialModP from a PolynomialModP where the two have different prime
p's.
"""
function poly_mod_p_invalid_subtract(; seed::Int = 0)
    Random.seed!(seed)
    # try and subtract two PolynomialModPs that have different p's
    p1 = PolynomialModP(rand(PolynomialSparse), 17)
    p2 = PolynomialModP(rand(PolynomialSparse), 5)
    try
        p1 - p2
    catch AssertionError
        println("poly_mod_p_invalid_subtract - PASSED")
        return
    end
    throw("PolynomialModP's that have different p's shouldn't be able to be subtracted" *
          "together")
end

"""
Test multiplying two PolynomialModPs that have different prime p's.
"""
function poly_mod_p_invalid_prod(; seed::Int = 0)
    Random.seed!(seed)
    # try and multiply two PolynomialModPs that have different p's
    p1 = PolynomialModP(rand(PolynomialSparse), 17)
    p2 = PolynomialModP(rand(PolynomialSparse), 5)
    try
        p1 * p2
    catch AssertionError
        println("poly_mod_p_invalid_prod - PASSED")
        return
    end
    throw("PolynomialModP's that have different p's shouldn't be able to be multipled " *
          "together")
end

"""
Test creating a PolynomialModP with a non-prime integer.
"""
function poly_mod_p_invalid_p(; seed::Int = 0)
    Random.seed!(seed)
    # try and give a composite number to a PolynomialModP constructor
    composite = rand(Int)
    while composite < 0 || isprime(composite)
        composite = rand(Int)
    end
    try
        p1 = PolynomialModP(rand(PolynomialSparse), composite)
    catch AssertionError
        println("poly_mod_p_invalid_p - PASSED")
        return
    end
    throw("PolynomialModP should not be able to take a composite number")
end

"""
Test that two PolynomialModPs with the same polynomial but different primes are not equal.
"""
function poly_mod_p_same_poly_diff_prime(; seed::Int = 0)
    Random.seed!(seed)
    poly_sparse = rand(PolynomialSparse)
    p1 = PolynomialModP(poly_sparse, 7)
    p2 = PolynomialModP(poly_sparse, 5)
    @assert p1 != p2
end
