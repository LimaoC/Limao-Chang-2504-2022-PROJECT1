#############################################################################
#############################################################################
#
# This file implements polynomial multiplication 
#                                                                               
#############################################################################
#############################################################################

"""
Chinese remainder theorem for two polynomials.
"""
function poly_crt(a::Polynomial, b::Polynomial, n::Integer, m::Integer)::PolynomialSparseBI
    @assert gcdx(n, m)[1] == 1
    c = PolynomialSparseBI()
    a_index = 1  # pointer to term of degree k in a
    b_index = 1  # pointer to term of degree k in b

    for k in 0:max(degree(a), degree(b))
        # find the term of degree k in a and b; if it doesn't exist, stop at the term that
        # would follow it
        # keep going if we aren't at the end of the array or we haven't reached the index
        while max(a_index, k) <= degree(a)  
            a.terms[a_index].degree >= k && break
            a_index += 1
        end
        while b_index <= degree(b) && k <= degree(b)
            b.terms[b_index].degree >= k && break
            b_index += 1
        end

        # get the coefficient of x^k in a and b
        ak = (k > degree(a) || iszero(a) || a.terms[a_index].degree != k) ?
             0 : a.terms[a_index].coeff
        bk = (k > degree(b) || iszero(b) || b.terms[b_index].degree != k) ?
             0 : b.terms[b_index].coeff
        ck = int_crt([ak, bk], [n, m])  # integer chinese remainder theorem
        if ck != 0
            c += Term(ck, k)
        end
    end

    return c
end

"""
Multiply two polynomials.
"""
function *(p1::P, p2::P)::Polynomial where {P<:Polynomial}
    p_out = P()
    for t in p1
        p_out += t * p2
    end
    return p_out
end
function old_mult(p1::PolynomialSparseBI, p2::PolynomialSparseBI)::PolynomialSparseBI
    # old multiplication function for benchmarking purposes
    p_out = PolynomialSparseBI()
    for t in p1
        p_out += t * p2
    end
    return p_out
end
function *(p1::PolynomialSparseBI, p2::PolynomialSparseBI)::PolynomialSparseBI
    (iszero(p1) || iszero(p2)) && return PolynomialSparseBI()

    height_p1 = maximum(map(term -> abs(term.coeff), p1.terms))
    height_p2 = maximum(map(term -> abs(term.coeff), p2.terms))

    # calculate bound for product of primes
    B = big(2) * height_p1 * height_p2 * min(length(p1.terms) + 1, length(p2.terms) + 1)
    p = 3
    M = big(p)  # product of primes
    c = (PolynomialModP(mod(p1, p), p) * PolynomialModP(mod(p2, p), p)).polynomial

    # repeatedly calculate poly_crt's until we have hit the bound B
    while M < B
        p = nextprime(p+1)
        c_prime = (PolynomialModP(mod(p1, p), p) * PolynomialModP(mod(p2, p), p)).polynomial
        c = poly_crt(c, c_prime, M, p)
        M *= p
    end

    return smod(c, M)
end
function *(p1::PolynomialModP, p2::PolynomialModP)
    @assert p1.prime == p2.prime
    return PolynomialModP(mod(p1.polynomial * p2.polynomial, p1.prime), p1.prime)
end

"""
Power of a polynomial.
"""
function ^(p::Polynomial, n::Int)::Polynomial
    n < 0 && error("No negative power")
    n == 0 && return one(p)

    out = one(p)
    squares = p
    
    # get truncated binary representation of exponent (i.e., most significant bit in string
    # is a 1)
    n_bin = reverse(bitstring(n)[findfirst('1', bitstring(n)):end])
    n_bin_length = length(n_bin)
    # iterate through in reverse order
    for (i, bit) in enumerate(n_bin)
        # square each iteration, and if bit is 1, multiply out by current value of squares
        if parse(Int, bit) == 1
            out *= squares
        end
        # don't need to square on the last iteration
        if i != n_bin_length
            squares *= squares
        end
    end
    return out
end
function ^(p::PolynomialModP, n::Int)::PolynomialModP
    n < 0 && error("No negative power")
    n == 0 && return PolynomialModP(one(PolynomialSparse), p.prime)

    out = one(PolynomialSparse)
    squares = p.polynomial
    
    # get truncated binary representation of exponent (i.e., most significant bit in string
    # is a 1)
    n_bin = reverse(bitstring(n)[findfirst('1', bitstring(n)):end])
    n_bin_length = length(n_bin)
    # iterate through in reverse order
    for (i, bit) in enumerate(n_bin)
        # square each iteration, and if bit is 1, multiply out by current value of squares
        if parse(Int, bit) == 1
            out = mod(out * squares, p.prime)
        end
        # don't need to square on the last iteration
        if i != n_bin_length
            squares = mod(squares * squares, p.prime)
        end
    end
    return PolynomialModP(out, p.prime)
end