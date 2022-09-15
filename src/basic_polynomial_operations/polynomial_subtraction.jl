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

"""
Subtract an integer from a polynomial.
"""
-(p::Polynomial, n::Int) = p + Term(-n, 0)
-(n::Int, p::Polynomial) = Term(n, 0) - p
