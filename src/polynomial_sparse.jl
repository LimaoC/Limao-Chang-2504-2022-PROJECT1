#############################################################################
#############################################################################
#
# This file defines the sparse polynomial type with several operations 
#                                                                               
#############################################################################
#############################################################################

####################################
# Polynomial type and construction #
####################################

"""
A PolynomialSparse type - designed to be for sparse polynomials with integer coefficients.
"""
struct PolynomialSparse <: Polynomial
    # A zero packed vector of terms
    # Only non-zero terms (i.e., terms with zero coefficient) are stored, and terms are
    # assumed to be in order with first term having degree 0, second degree 1, and so forth
    # until the degree of the polynomial. The last term is assumed to be non-zero except
    # for the zero polynomial where the vector is of length 1.
    terms::Vector{Term{Int64}}

    # Inner constructor of 0 polynomial
    PolynomialSparse() = new([])

    # Inner constructor of polynomial based on arbitrary list of terms
    function PolynomialSparse(vt::Vector{Term{Int64}})
        # Filter the vector so that there is not more than a single zero term
        vt = filter((t) -> !iszero(t), vt)
        if isempty(vt)
            vt = [zero(Term)]
        end

        terms = filter((t) -> t.coeff != 0, vt)  # filter out zero terms
        return new(terms)
    end
end

"""
Construct a (sparse) polynomial with a single term.
"""
PolynomialSparse(t::Term) = PolynomialSparse([t])

################################
# Pushing and popping of terms #
################################

"""
Push a new term into the polynomial.
"""
# Note that ideally this would throw and error if pushing another term of degree that is
# already in the polynomial
function push!(p::PolynomialSparse, t::Term)
    if t.degree > degree(p)
        push!(p.terms, t)
    else
        # find where we should store this term
        index = findfirst(term -> term.degree >= t.degree, p.terms)
        if index == nothing
            # insertion term is smaller than any other term we have; insert at start
            insert!(p.terms, 1, t)
        elseif t.degree < p.terms[index].degree
            # term of this degree doesn't yet exist so insert a new term
            insert!(p.terms, index, t)
        else
            # term of this degree already exists
            p.terms[index] = t
        end
    end
    return p
end

"""
Pop the leading term out of the polynomial.
"""
function pop!(p::PolynomialSparse)::Term
    popped_term = pop!(p.terms)  # last element popped is leading coefficient

    if isempty(p.terms)
        push!(p.terms, zero(Term))
    end

    return popped_term
end
