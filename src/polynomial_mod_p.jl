#############################################################################
#############################################################################
#
# This file defines the PolynomialModP type with several operations 
#                                                                               
#############################################################################
#############################################################################

########################################
# PolynomialModP type and construction #
########################################

"""
A PolynomialModP type - holds a Polynomial and a prime number.
"""
struct PolynomialModP <: Polynomial
    polynomial::PolynomialSparse
    prime::Integer

    # Inner constructor of PolynomialModP
    function PolynomialModP(p::PolynomialSparse, prime::Integer)
        return new(p, prime)
    end
end

###########
# Display #
###########
function show(io::IO, p::PolynomialModP)
    show(io, p.polynomial)
    print(io, " (mod $(p.prime))")
end

################################
# Pushing and popping of terms #
################################

"""
Push a new term into the polynomial.
"""
push!(p::PolynomialModP, t::Term) = push!(p.polynomial, t)

"""
Pop the leading term out of the polynomial.
"""
pop!(p::PolynomialModP)::Term = pop!(p.polynomial)

"""
Power of a polynomial mod prime.
"""
pow_mod(p::PolynomialSparse, n::Int, prime::Int) = PolynomialModP(p, prime)^n