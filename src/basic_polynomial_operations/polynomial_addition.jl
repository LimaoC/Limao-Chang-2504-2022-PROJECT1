#############################################################################
#############################################################################
#
# This file implements polynomial addition 
#                                                                               
#############################################################################
#############################################################################

"""
Add a polynomial and a term.
"""
function +(p::PolynomialDense, t::Term)
    p = deepcopy(p)
    if t.degree > degree(p)
        push!(p, t)
    else
        if !iszero(p.terms[t.degree+1])  # + 1 is due to indexing
            p.terms[t.degree+1] += t
        else
            p.terms[t.degree+1] = t
        end
    end

    return trim!(p)
end
function +(p::PolynomialSparse, t::Term)
    p = deepcopy(p)
    if t.degree > degree(p)
        push!(p, t)
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
            # term of this degree already exists so add to existing term
            p.terms[index] += t
        end
    end
    return trim!(p)
end
function +(p::PolynomialSparseBI, t::Term)
    t = Term(big(t.coeff), t.degree)  # convert to Term{BigInt} before addition
    p = deepcopy(p)
    if t.degree > degree(p)
        push!(p, t)
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
            # term of this degree already exists so add to existing term
            p.terms[index] += t
        end
    end
    return trim!(p)
end
+(t::Term, p::Polynomial) = p + t

"""
Add two polynomials.
"""
function +(p1::Polynomial, p2::Polynomial)::Polynomial
    p = deepcopy(p1)
    for t in p2
        p += t
    end
    return p
end

"""
Add a polynomial and an integer.
"""
+(p::Polynomial, n::Int) = p + Term(n, 0)
+(n::Int, p::Polynomial) = p + Term(n, 0)