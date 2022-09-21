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
    a_index = 1
    b_index = 1
    x = x_poly(PolynomialSparseBI)

    for k in 0:max(degree(a), degree(b))
        # find the term of degree k in a and b; if it doesn't exist, stop at the term that
        # would follow it
        while a_index <= degree(a) && k <= degree(a)
            if a.terms[a_index].degree >= k
                break
            else
                a_index += 1
            end
        end
        while b_index <= degree(b) && k <= degree(b)
            if b.terms[b_index].degree >= k
                break
            else
                b_index += 1
            end
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
        new_summand = (t * p2)
        p_out = p_out + new_summand
    end
    return p_out
end
function *(p1::PolynomialSparseBI, p2::PolynomialSparseBI)::PolynomialSparseBI
    if iszero(p1) || iszero(p2)
        return PolynomialSparseBI()
    end

    height_p1 = maximum(map(term -> abs(term.coeff), p1.terms))
    height_p2 = maximum(map(term -> abs(term.coeff), p2.terms))

    B = 2 * height_p1 * height_p2 * min(length(p1.terms) + 1, length(p2.terms) + 1)
    p = 3
    M = big(p)
    c = PolynomialModP(mod(p1, p), p) * PolynomialModP(mod(p2, p), p)

    while M < B
        p = nextprime(p+1)
        c_prime = PolynomialModP(mod(p1, p), p) * PolynomialModP(mod(p2, p), p)
        c = poly_crt(c, c_prime, M, p)
        M *= p
    end

    return smod(c, M)
end
function *(p1::PolynomialModP, p2::PolynomialModP)
    @assert p1.prime == p2.prime
    return mod(p1.polynomial * p2.polynomial, p1.prime)
end

"""
Power of a polynomial.
"""
function ^(p::Polynomial, n::Int)
    n < 0 && error("No negative power")
    out = one(p)
    for _ in 1:n
        out *= p
    end
    return out
end
^(p::PolynomialModP, n::Int) = begin
    return mod(p.polynomial^n, p.prime)
end