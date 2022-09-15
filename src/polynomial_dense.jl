#############################################################################
#############################################################################
#
# This file defines the dense polynomial type with several operations 
#                                                                               
#############################################################################
#############################################################################

####################################
# Polynomial type and construction #
####################################

"""
A PolynomialDense type - designed to be for dense polynomials with integer coefficients.
Stores terms densely - i.e. all terms up to the highest degree are stored.
"""
struct PolynomialDense <: Polynomial
    # A zero packed vector of terms
    # Terms are assumed to be in order with first term having degree 0, second degree 1,
    # and so forth until the degree of the polynomial. The leading term (i.e. last) is
    # assumed to be non-zero except for the zero polynomial where the vector is of length 1.
    # Note: at positions where the coefficient is 0, the power of the term is also 0 (this
    # is how the Term type is designed)
    terms::Vector{Term{Int64}}

    # Inner constructor of 0 polynomial
    PolynomialDense() = new([zero(Term)])

    # Inner constructor of polynomial based on arbitrary list of terms
    function PolynomialDense(vt::Vector{Term{Int64}})
        # Filter the vector so that there is not more than a single zero term
        vt = filter((t) -> !iszero(t), vt)
        if isempty(vt)
            vt = [zero(Term)]
        end

        max_degree = maximum((t) -> t.degree, vt)
        terms = [zero(Term) for i in 0:max_degree]  # set all terms to zero to start with

        # update based on input terms
        for t in vt
            terms[t.degree+1] = t  # + 1 accoutns for 1-indexing
        end

        return new(terms)
    end
end

"""
Construct a (dense) polynomial with a single term.
"""
PolynomialDense(t::Term) = PolynomialDense([t])

################################
# Pushing and popping of terms #
################################

"""
Push a new term into the polynomial.
"""
# Note that ideally this would throw and error if pushing another term of degree that is
# already in the polynomial
function push!(p::PolynomialDense, t::Term)
    if t.degree <= degree(p)
        p.terms[t.degree+1] = t
    else
        append!(p.terms, zeros(Term, t.degree - degree(p) - 1))
        push!(p.terms, t)
    end
    return p
end

"""
Pop the leading term out of the polynomial. When polynomial is 0, keep popping out 0.
"""
function pop!(p::PolynomialDense)::Term
    popped_term = pop!(p.terms)  # last element popped is leading coefficient

    while !isempty(p.terms) && iszero(last(p.terms))
        pop!(p.terms)
    end

    if isempty(p.terms)
        push!(p.terms, zero(Term))
    end

    return popped_term
end