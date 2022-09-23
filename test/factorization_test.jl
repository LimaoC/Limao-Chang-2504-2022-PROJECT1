#############################################################################
#############################################################################
#
# This file contains units tests for polynomial factorization
#                                                                               
#############################################################################
#############################################################################

"""
Executes all polynomial factorization tests in this file.
"""
function factorization_tests()
    factor_test_poly_dense()
    factor_test_poly_sparse()
    factor_test_poly_sparse_bi()
end

"""
Test factorization of dense polynomials.
"""
function factor_test_poly_dense(;N::Int = 10, seed::Int = 0,
                                primes::Vector{Int} = [5,17,19])
    Random.seed!(seed)
    for prime in primes
        print("\ndoing prime = $prime \t")
        for _ in 1:N
            print(".")
            p = rand(PolynomialDense)
            factorization = factor(p, prime)
            pr = mod(expand_factorization(factorization),prime)
            @assert mod(p-pr,prime) == 0 
        end
    end

    println("\nfactor_test_poly_dense - PASSED")
end

"""
Test factorization of sparse polynomials.
"""
function factor_test_poly_sparse(;N::Int = 10, seed::Int = 0, primes::Vector{Int} = [5,17,19])
    Random.seed!(seed)
    for prime in primes
        print("\ndoing prime = $prime \t")
        for _ in 1:N
            print(".")
            p = rand(PolynomialSparse)
            factorization = factor(p, prime)
            pr = mod(expand_factorization(factorization),prime)
            @assert mod(p-pr,prime) == 0 
        end
    end

    println("\nfactor_test_poly_sparse - PASSED")
end

"""
Test factorization of sparse big int polynomials.
"""
function factor_test_poly_sparse_bi(;N::Int = 10, seed::Int = 0, primes::Vector{Int} = [5,17,19])
    Random.seed!(seed)
    for prime in primes
        print("\ndoing prime = $prime \t")
        for _ in 1:N
            print(".")
            p = rand(PolynomialSparse)
            factorization = factor(p, prime)
            pr = mod(expand_factorization(factorization),prime)
            @assert mod(p-pr,prime) == 0 
        end
    end

    println("\nfactor_test_poly_sparse_bi - PASSED")
end
