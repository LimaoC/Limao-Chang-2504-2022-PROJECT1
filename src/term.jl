#############################################################################
#############################################################################
#
# This file defines the Term type with several operations 
#                                                                               
#############################################################################
#############################################################################

##############################
# Term type and construction #
##############################

"""
A term.
"""
struct Term{T<:Integer}  # structs are immutable by default
    coeff::T
    degree::Int64
    function Term(coeff::T, degree::Int64) where {T<:Integer}
        degree < 0 && error("Degree must be non-negative")
        coeff != 0 ? new{T}(coeff, degree) : new{T}(coeff, 0)
    end
end

"""
Creates the zero term. If no parametric type is specified, the type defaults to Int64.
"""
zero(::Type{Term})::Term = Term(0, 0)
(zero(::Type{Term{T}})::Term) where {T<:Integer} = Term(T(0), 0)

"""
Creates the unit term. If no parametric type is specified, the type defaults to Int64.
"""
one(::Type{Term})::Term = Term(1, 0)
(one(::Type{Term{T}})::Term) where {T<:Integer} = Term(T(1), 0)

###########
# Display #
###########

"""
Show a term.
"""
show(io::IO, t::Term) = begin
    # omit degree for terms with a degree of 1, omit x for constant terms
    x = t.degree > 1 ? "x^$(t.degree)" : (t.degree == 1 ? "x" : "")
    # omit coefficient if it is ±1, except if it is a constant term
    coefficient = ((t.degree > 0 && abs(t.coeff) != 1) || t.degree == 0) ?
                  t.coeff : (t.coeff == 1 ? "" : "-")
    print(io, "$coefficient$x")
end

########################
# Queries about a term #
########################

"""
Check if a term is 0.
"""
iszero(t::Term)::Bool = iszero(t.coeff)

"""
Compare two terms.
"""
isless(t1::Term, t2::Term)::Bool = begin
    t1.degree == t2.degree ? (t1.coeff < t2.coeff) : (t1.degree < t2.degree)
end

"""
Evaluate a term at a point x.
"""
evaluate(t::Term, x::T) where {T<:Number} = t.coeff * x^t.degree

###########################
# Queries about two terms #
###########################

"""
Check if two terms are equal.
"""
==(t1::Term, t2::Term)::Bool = t1.coeff == t2.coeff && t1.degree == t2.degree

##########################
# Operations with a term #
##########################

"""
Add two terms of the same degree.
"""
function +(t1::Term, t2::Term)::Term
    @assert t1.degree == t2.degree
    Term(t1.coeff + t2.coeff, t1.degree)
end

"""
Negate a term.
"""
-(t::Term,) = Term(-t.coeff, t.degree)

"""
Subtract two terms with the same degree.
"""
function -(t1::Term, t2::Term)::Term
    @assert t1.degree == t2.degree
    return t1 + (-t2)
end

"""
Multiply two terms.
"""
*(t1::Term, t2::Term)::Term = Term(t1.coeff * t2.coeff, t1.degree + t2.degree)

"""
Compute the mod of a term with an integer.
"""
mod(t::Term, p::Integer) = Term(mod(t.coeff, p), t.degree)

"""
Symmetric mod.
"""
smod(t::Term, p::Integer)::Term = Term(smod(t.coeff, p), t.degree)

"""
Compute the derivative of a term.
"""
derivative(t::Term) = Term(t.coeff * t.degree, max(t.degree - 1, 0))

"""
Divide two terms. Returns a function of an integer.
"""
function ÷(t1::Term, t2::Term)
    @assert t1.degree ≥ t2.degree
    f(p::Int)::Term = Term(mod((t1.coeff * int_inverse_mod(t2.coeff, p)), p),
                           t1.degree - t2.degree)
end

"""
Integer divide a term by an integer.
"""
÷(t::Term, n::Int) = t ÷ Term(n, 0)
