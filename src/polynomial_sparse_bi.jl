#############################################################################
#############################################################################
#
# This file defines the sparse polynomial (big int) type with several
# operations 
#                                                                               
#############################################################################
#############################################################################

####################################
# Polynomial type and construction #
####################################

"""
A PolynomialSparseBI type - designed to be for sparse polynomials with big integer
coefficients.
"""
struct PolynomialSparseBI <: Polynomial
    # A zero packed vector of terms
    # Only non-zero terms (i.e., terms with zero coefficient) are stored, and terms are
    # assumed to be in order with first term having degree 0, second degree 1, and so forth
    # until the degree of the polynomial. The last term is assumed to be non-zero except
    # for the zero polynomial where the vector is of length 1.
    terms::Vector{Term{BigInt}}

    # Inner constructor of 0 polynomial
    PolynomialSparseBI() = new([])

    # Inner constructor of polynomial based on arbitrary list of terms
    function PolynomialSparseBI(vt::Vector{Term{BigInt}})
        # Filter the vector so that there is not more than a single zero term
        vt = filter((t) -> !iszero(t), vt)
        if isempty(vt)
            vt = [zero(Term)]
        end

        terms = filter((t) -> t.coeff != 0, vt)  # filter out zero terms
        return new(terms)
    end

    function PolynomialSparseBI(vt::Vector{Term{T}}) where {T<:Integer}
        vt = map((t) -> Term(big(t.coeff), t.degree), vt)
        return PolynomialSparseBI(vt)
    end
end

"""
Construct a (sparse, big integer) polynomial with a single term.
"""
PolynomialSparseBI(t::Term) = PolynomialSparseBI([t])

################################
# Pushing and popping of terms #
################################

"""
Push a new term into the polynomial.
"""
# Note that ideally this would throw and error if pushing another term of degree that is
# already in the polynomial
function push!(p::PolynomialSparseBI, t::Term)
    t = Term(big(t.coeff), t.degree)  # convert to Term{BigInt} before pushing
    if t.degree > degree(p)
        push!(p.terms, t)
    else
        index = findfirst(term -> term.degree >= t.degree, p.terms)
        if index != nothing
            p.terms[index] = t
        else
            insert!(p.terms, 1, t)
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

##################################################################
# Operations with two objects where at least one is a polynomial #
##################################################################

"""
Take the mod of a sparse big int polynomial with an integer and convert to machine int.
"""
function mod(f::PolynomialSparseBI, p::Integer)::PolynomialSparse
    f_out = PolynomialSparse()
    for i in 1:length(f.terms)
        term = mod(f.terms[i], p)
        !iszero(term) && push!(f_out, term)
    end
    return trim!(f_out)

    # p_out = Polynomial()
    # for t in f
    #     new_term = mod(t, p)
    #     @show new_term
    #     push!(p_out, new_term)
    # end
    # return p_out
end

"""
Take the symmetric mod of a sparse big int polynomial with an integer.
"""
function smod(f::PolynomialSparseBI, p::Integer)::PolynomialSparseBI
    f_out = PolynomialSparseBI()
    for i in 1:length(f.terms)
        push!(f_out, smod(f.terms[i], p))
    end
    return trim!(f_out)
end
