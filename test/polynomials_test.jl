#############################################################################
#############################################################################
#
# This file contains units tests for polynomial operations
#                                                                               
#############################################################################
#############################################################################


"""
Test product of polynomials.
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
Test derivative of polynomials (as well as product).
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

"""
Test division of polynomials modulo p.
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

"""
Test the extended euclid algorithm for polynomials modulo p.
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