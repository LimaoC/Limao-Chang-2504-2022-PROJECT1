#############################################################################
#############################################################################
#
# This file implements the polynomial subtraction that isn't covered by
# polynomial addition
#                                                                               
#############################################################################
#############################################################################
"""
Subtract a polynomial from a term.
"""
-(t::Term, p::Polynomial) = t + (-p)
-(t::Term, p::Polynomial) = t + (-p.polynomial)

"""
Subtract an integer from a polynomial.
"""
-(p::Polynomial, n::Int) = p + Term(-n, 0)
-(n::Int, p::Polynomial) = Term(n, 0) - p
-(p::PolynomialModP, n::Int) = p.polynomial + Term(-n, 0)
-(p::PolynomialModP, n::Int) = Term(n, 0) - p.polynomial
