#############################################################################
#############################################################################
#
# This file defines the abstract polynomial type with several operations 
#                                                                               
#############################################################################
#############################################################################

####################################
# Polynomial type and construction #
####################################

"""
A Polynomial type - designed to be for polynomials with integer coefficients.
"""
abstract type Polynomial end

"""
This function maintains the invariant of the Polynomial type so that there are no zero
terms beyond the highest non-zero term.
"""
function trim!(p::Polynomial)::Polynomial
    i = length(p.terms)
    while i > 1
        if iszero(p.terms[i])
            pop!(p.terms)
        else
            break
        end
        i -= 1
    end
    return p
end

"""
Construct a polynomial of the form x^p-x.
"""
(cyclotonic_polynomial(::Type{P}, p::Int)::P) where {P<:Polynomial} = begin
    P([Term(1, p), Term(-1, 0)])
end

"""
Construct a polynomial of the form x-n.
"""
linear_monic_polynomial(::Type{P}, n::Int) where {P<:Polynomial} = begin
    P([Term(1, 1), Term(-n, 0)])
end

"""
Construct a polynomial of the form x.
"""
x_poly(::Type{P}) where {P<:Polynomial} = P(Term(1, 1))

"""
Creates the zero polynomial.
"""
# Note: in a parametric function, brackets need to wrap the function declaration if we want
# to specify a return type.
# Reference: https://discourse.julialang.org/t/how-to-declare-parametric-with-return-type-of-a-function/40288
(zero(::Type{P})::P) where {P<:Polynomial} = P()

"""
Creates the unit polynomial.
"""
(one(::Type{P})::P) where {P<:Polynomial} = P(one(Term))
one(p::Polynomial) = one(typeof(p))

"""
Generates a random polynomial.
"""
function rand(::Type{P};
    degree::Int=-1,
    terms::Int=-1,
    max_coeff::Int=100,
    mean_degree::Float64=5.0,
    prob_term::Float64=0.7,
    monic=false,
    condition=(p) -> true)::P where {P<:Polynomial}

    while true
        _degree = degree == -1 ? rand(Poisson(mean_degree)) : degree
        _terms = terms == -1 ? rand(Binomial(_degree, prob_term)) : terms
        degrees = vcat(sort(sample(0:_degree-1, _terms, replace=false)), _degree)
        coeffs = rand(1:max_coeff, _terms + 1)
        monic && (coeffs[end] = 1)
        p = P([Term(coeffs[i], degrees[i]) for i in 1:length(degrees)])
        condition(p) && return p
    end
end

###########
# Display #
###########

"""
Show a polynomial.
"""
function show(io::IO, p::Polynomial)
    if iszero(p)
        print(io, "0")
    else
        n = length(p.terms)
        # print terms in descending degree order
        for (i, t) in enumerate(reverse(p.terms))
            if !iszero(t)
                # if first term, omit ± signs
                # if negative coefficient, detach negative sign from term 
                print(io, i == 1 ? t : (t.coeff < 0 ? " - $(string(t)[2:end])" : " + $t"))
            end
        end
    end
end

##############################################
# Iteration over the terms of the polynomial #
##############################################

"""
Allows to do iteration over the non-zero terms of the polynomial. This implements the
iteration interface.
"""
iterate(p::Polynomial, state=1) = iterate(p.terms, state)

##############################
# Queries about a polynomial #
##############################

"""
The number of terms of the polynomial.
"""
length(p::Polynomial) = length(p.terms)

"""
The leading term of the polynomial.
"""
leading(p::Polynomial)::Term = isempty(p.terms) ? zero(Term) : last(p.terms)

"""
Returns the coefficients of the polynomial.
"""
coeffs(p::Polynomial)::Vector{Int} = [t.coeff for t in p]

"""
The degree of the polynomial.
"""
degree(p::Polynomial)::Int = leading(p).degree

"""
The content of the polynomial is the GCD of its coefficients.
"""
content(p::Polynomial)::Int = euclid_alg(coeffs(p))

"""
Evaluate the polynomial at a point `x`.
"""
evaluate(f::Polynomial, x::T) where {T<:Number} = sum(evaluate(t, x) for t in f)

"""
Check if the polynomial is zero.
"""
iszero(p::Polynomial)::Bool = p.terms == [Term(0, 0)] || isempty(p.terms)

#################################################################
# Transformation of the polynomial to create another polynomial #
#################################################################

"""
The negative of a polynomial.
"""
-(p::P) where {P<:Polynomial} = P(map((pt) -> -pt, p.terms))

"""
Create a new polynomial which is the derivative of the polynomial.
"""
function derivative(p::P)::P where {P<:Polynomial}
    der_p = P()
    for term in p
        der_term = derivative(term)
        !iszero(der_term) && push!(der_p, der_term)
    end
    return trim!(der_p)
end

"""
The prim part (multiply a polynomial by the inverse of its content).
"""
prim_part(p::Polynomial) = p ÷ content(p)

"""
A square free polynomial.
"""
square_free(p::Polynomial, prime::Int)::Polynomial = begin
    (p ÷ gcd(p, derivative(p), prime))(prime)
end

#################################
# Queries about two polynomials #
#################################

"""
Check if two polynomials are the same.
"""
(==(p1::P, p2::P)::Bool) where {P<:Polynomial} = p1.terms == p2.terms

"""
Check if a polynomial is equal to 0.
"""
# Note that in principle there is a problem here. E.g The polynomial 3 will return true to
# equalling the integer 2.
==(p::Polynomial, n::T) where {T<:Real} = iszero(p) == iszero(n)

##################################################################
# Operations with two objects where at least one is a polynomial #
##################################################################

"""
Subtraction of two polynomials.
"""
-(p1::Polynomial, p2::Polynomial)::Polynomial = p1 + (-p2)

"""
Multiplication of polynomial and term.
"""
function *(t::Term, p1::P)::P where {P<:Polynomial}
    iszero(t) ? P() : P(map((pt) -> t * pt, p1.terms))
end
*(p1::Polynomial, t::Term)::Polynomial = t * p1

"""
Multiplication of polynomial and an integer.
"""
*(n::Integer, p::Polynomial)::Polynomial = p * Term(n, 0)
*(p::Polynomial, n::Integer)::Polynomial = n * p

"""
Integer division of a polynomial by an integer.

Warning this may not make sense if n does not divide all the coefficients of p.
"""
function ÷(p::P, n::Int) where {P<:Polynomial}
    (prime) -> P(map((pt) -> ((pt ÷ n)(prime)), p.terms))
end

"""
Take the mod of a polynomial with an integer.
"""
function mod(f::P, p::Integer)::P where {P<:Polynomial}
    f_out = P()
    for i in 1:length(f.terms)
        term = mod(f.terms[i], p)
        !iszero(term) && push!(f_out, term)
    end
    return trim!(f_out)
end

"""
Power of a polynomial mod prime.
"""
function pow_mod(p::Polynomial, n::Int, prime::Int)::Polynomial
    n < 0 && error("No negative power")
    out = one(p)
    for _ in 1:n
        out *= p
        out = mod(out, prime)
    end
    return out
end
