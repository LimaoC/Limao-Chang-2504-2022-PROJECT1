#############################################################################
#############################################################################
#
# This file contains unit tests for polynomial pow_mod and ^ operations
#                                                                               
#############################################################################
#############################################################################

function polynomial_power_tests()
    poly_sparse_bi_pow_test()
    poly_mod_p_pow_test()
    poly_sparse_bi_pow_mod_test()
    poly_mod_p_pow_mod_test()
end

"""
Test whether ^() raises the power of a PolynomialSparseBI correctly.
"""
function poly_sparse_bi_pow_test(; seed::Int = 0, N::Int = 10, n::Int = 10)
    Random.seed!(seed)
    for _ in 1:N
        p = rand(PolynomialSparseBI)
        p_out = deepcopy(p)
        for _ in 1:n-1
            p_out *= p
        end
        @assert p^n == p_out
    end
    println("poly_sparse_bi_pow_test - PASSED")
end

"""
Test whether ^() raises the power of a PolynomialModP correctly.
"""
function poly_mod_p_pow_test(; seed::Int = 0, N::Int = 10, n::Int = 10)
    Random.seed!(seed)
    for _ in 1:N
        p = rand(PolynomialModP)
        p_out = deepcopy(p)
        for _ in 1:n-1
            p_out *= p
        end
        @assert p^n == p_out
    end
    println("poly_mod_p_pow_test - PASSED")
end

"""
Test whether pow_mod() raises the power of a PolynomialSparseBI correctly.
"""
function poly_sparse_bi_pow_mod_test(; seed::Int = 0, N::Int = 10, n::Int = 10,
                                     prime_num::Int = 17)
    Random.seed!(seed)
    for _ in 1:N
        p = rand(PolynomialSparseBI)
        p_out = deepcopy(p)
        for _ in 1:n-1
            p_out = mod(p_out * p, prime_num)
        end
        @assert pow_mod(p, n, prime_num) == p_out
    end
    println("poly_sparse_bi_pow_mod_test - PASSED")
end

"""
Test whether pow_mod() raises the power of a PolynomialModP correctly.
"""
function poly_mod_p_pow_mod_test(; seed::Int = 0, N::Int = 10, n::Int = 10,
                                     prime_num::Int = 17)
    Random.seed!(seed)
    for _ in 1:N
        p = rand(PolynomialModP, p=prime_num)
        p_out = deepcopy(p)
        for _ in 1:n-1
            p_out *= p
        end
        @assert pow_mod(p, n) == p_out
    end
    println("poly_mod_p_pow_mod_test - PASSED")
end
