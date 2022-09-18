#############################################################################
#############################################################################
#
# This file implements polynomial multiplication 
#                                                                               
#############################################################################
#############################################################################

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
*(p1::PolynomialModP, p2::PolynomialModP) = begin
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