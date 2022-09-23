#############################################################################
#############################################################################
#
# This file contains units tests for polynomial division operations
#                                                                               
#############################################################################
#############################################################################

"""
Executes all polynomial division tests in this file
"""
function polynomial_division_tests()
    division_test_poly_dense()
    division_test_poly_sparse()
    division_test_poly_sparse_bi()
end

"""
Test division of dense polynomials modulo p.
"""
function division_test_poly_dense(; prime::Int=101, N::Int=10^4, seed::Int=0)
    Random.seed!(seed)
    for _ in 1:N
        p1 = rand(PolynomialDense)
        p2 = rand(PolynomialDense)
        p_prod = p1 * p2
        q, r = PolynomialDense(), PolynomialDense()
        try
            q, r = divide(p_prod, p2)(prime)
            if (q, r) == (nothing, nothing)
                println("Unlucky prime: $p1 is reduced to $(p1 % prime) modulo $prime")
                continue
            end
        catch e
            if typeof(e) == DivideError
                @assert mod(p2, prime) == 0
            else
                throw(e)
            end
        end
        @assert iszero(mod(q * p2 + r - p_prod, prime))
    end
    println("division_test_poly_dense - PASSED")
end

"""
Test division of sparse polynomials modulo p.
"""
function division_test_poly_sparse(; prime::Int=101, N::Int=10^4, seed::Int=0)
    Random.seed!(seed)
    for _ in 1:N
        p1 = rand(PolynomialSparse)
        p2 = rand(PolynomialSparse)
        p_prod = p1 * p2
        q, r = PolynomialSparse(), PolynomialSparse()
        try
            q, r = divide(p_prod, p2)(prime)
            if (q, r) == (nothing, nothing)
                println("Unlucky prime: $p1 is reduced to $(p1 % prime) modulo $prime")
                continue
            end
        catch e
            if typeof(e) == DivideError
                @assert mod(p2, prime) == 0
            else
                throw(e)
            end
        end
        @assert iszero(mod(q * p2 + r - p_prod, prime))
    end
    println("division_test_poly_sparse - PASSED")
end

"""
Test division of sparse big int polynomials modulo p.
"""
function division_test_poly_sparse_bi(; prime::Int=101, N::Int=10^4, seed::Int=0)
    Random.seed!(seed)
    for _ in 1:N
        p1 = rand(PolynomialSparseBI)
        p2 = rand(PolynomialSparseBI)
        p_prod = p1 * p2
        q, r = PolynomialSparseBI(), PolynomialSparseBI()
        try
            q, r = divide(p_prod, p2)(prime)
            if (q, r) == (nothing, nothing)
                println("Unlucky prime: $p1 is reduced to $(p1 % prime) modulo $prime")
                continue
            end
        catch e
            if typeof(e) == DivideError
                @assert mod(p2, prime) == 0
            else
                throw(e)
            end
        end
        @assert iszero(mod(q * p2 + r - p_prod, prime))
    end
    println("division_test_poly_sparse_bi - PASSED")
end
