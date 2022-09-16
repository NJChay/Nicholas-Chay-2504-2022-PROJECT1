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
struct BTerm  #structs are immutable by default
    coeff::BigInt
    degree::Int
    function BTerm(coeff::BigInt, degree::Int)
        degree < 0 && error("Degree must be non-negative")
        coeff != 0 ? new(coeff,degree) : new(coeff,0)
    end
end

function coeff(t::BTerm)
    return t.coeff
end

function degree(t::BTerm)
    return t.degree
end

"""
Creates the zero term.
"""
zero(::Type{BTerm})::BTerm = BTerm(BigInt(0),0)

"""
Creates the unit term.
"""
one(::Type{BTerm})::BTerm = BTerm(BigInt(1),0)

###########
# Display #
###########

"""
Show a term.
"""
function show(io::IO, t::BTerm) 
    if t.degree!=0
         t.degree!=1 ? print(io, "$(t.coeff)⋅x^$(t.degree)") : print(io, "$(t.coeff)⋅x") 
    else print(io, "$(t.coeff)")#\cdot + [TAB]
    end
end

########################
# Queries about a term #
########################

"""
Check if a term is 0.
"""
iszero(t::BTerm)::Bool = iszero(t.coeff)

"""
Compare two terms.
"""
isless(t1::BTerm,t2::BTerm)::Bool =  t1.degree == t2.degree ? (t1.coeff < t2.coeff) : (t1.degree < t2.degree)  

"""
Evaluate a term at a point x.
"""
evaluate(t::BTerm, x::T) where T <: Number = t.coeff * x^t.degree

##########################
# Operations with a term #
##########################

"""
Add two terms of the same degree.
"""
function +(t1::BTerm,t2::BTerm)::BTerm
    @assert t1.degree == t2.degree
    BTerm(t1.coeff + t2.coeff, t1.degree)
end

"""
Negate a term.
"""
-(t::BTerm,) = BTerm(-t.coeff,t.degree)  

"""
Subtract two terms with the same degree.
"""
-(t1::BTerm, t2::BTerm)::BTerm = t1 + (-t2) 

"""
Multiply two terms.
"""
*(t1::BTerm, t2::BTerm)::BTerm = BTerm(t1.coeff * t2.coeff, t1.degree + t2.degree)

"""
Compute the symmetric mod of a term with an integer.
"""
mod(t::BTerm, p::Int) = BTerm(mod(t.coeff,p), t.degree)

"""
Compute the derivative of a term.
"""
derivative(t::BTerm) = BTerm(t.coeff*t.degree,max(t.degree-1,0))

"""
Divide two terms. Returns a function of an integer.
"""
function ÷(t1::BTerm,t2::BTerm) #\div + [TAB]
    @assert t1.degree ≥ t2.degree
    f(p::Int)::BTerm = BTerm(mod((t1.coeff * int_inverse_mod(t2.coeff, p)), p), t1.degree - t2.degree)
end

==(t1::BTerm, t2::BTerm) = t1.coeff==t2.coeff && t1.degree==t2.degree ? true : false

"""
Integer divide a term by an integer.
"""
÷(t::BTerm, n::Int) = t ÷ BTerm(BigInt(n),0)
÷(t::BTerm, n::BigInt) = t ÷ BTerm(n,0)
