#############################################################################
#############################################################################
#
# This file implements polynomial division 
#                                                                               
#############################################################################
#############################################################################

"""  Modular algorithm.
f divide by g

f = q*g + r

p is a prime
"""
function divide(num::PolynomialSparse, den::PolynomialSparse)
    function division_function(p::Int)
        f, g = mod(num,p), mod(den,p)
        degree(f) < degree(num) && return nothing 
        iszero(g) && throw(DivideError())
        q = PolynomialSparse()
        prev_degree = degree(f)
        while degree(f) ≥ degree(g) 
            h = PolynomialSparse( (leading(f) ÷ leading(g))(p) )  #syzergy 
            f = mod((f - h*g), p)
            q = mod((q + h), p)  
            prev_degree == degree(f) && break
            prev_degree = degree(f)
        end
        @assert iszero( mod((num  - (q*g + f)),p))
        return q, f
    end
    return division_function
end

function divide(num::PolynomialSparseBI, den::PolynomialSparseBI)
    function division_function(p::Int)
        f, g = mod(num,p), mod(den,p)
        degree(f) < degree(num) && return nothing 
        iszero(g) && throw(DivideError())
        q = PolynomialSparseBI()
        prev_degree = degree(f)
        while degree(f) ≥ degree(g) 
            h = PolynomialSparseBI( (leading(f) ÷ leading(g))(p) )  #syzergy 
            f = mod((f - h*g), p)
            q = mod((q + h), p)  
            prev_degree == degree(f) && break
            prev_degree = degree(f)
        end
        @assert iszero( mod((num  - (q*g + f)),p))
        return q, f
    end
    return division_function
end

function divide(num::PolynomialDense, den::PolynomialDense)
    function division_function(p::Int)
        f, g = mod(num,p), mod(den,p)
        degree(f) < degree(num) && return nothing 
        iszero(g) && throw(DivideError())
        q = PolynomialDense()
        prev_degree = degree(f)
        while degree(f) ≥ degree(g) 
            h = PolynomialDense( (leading(f) ÷ leading(g))(p) )  #syzergy 
            f = mod((f - h*g), p)
            q = mod((q + h), p)  
            prev_degree == degree(f) && break
            prev_degree = degree(f)
        end
        @assert iszero( mod((num  - (q*g + f)),p))
        return q, f
    end
    return division_function
end

"""
The quotient from polynomial division. Returns a function of an integer.
"""
÷(num::PolynomialSparse, den::PolynomialSparse)  = (p::Int) -> first(divide(num,den)(p))
÷(num::PolynomialDense, den::PolynomialDense)  = (p::Int) -> first(divide(num,den)(p))
÷(num::PolynomialSparseBI, den::PolynomialSparseBI)  = (p::Int) -> first(divide(num,den)(p))
"""
The remainder from polynomial division. Returns a function of an integer.
"""
rem(num::PolynomialSparse, den::PolynomialSparse)  = (p::Int) -> last(divide(num,den)(p))
rem(num::PolynomialDense, den::PolynomialDense)  = (p::Int) -> last(divide(num,den)(p))
rem(num::PolynomialSparseBI, den::PolynomialSparseBI)  = (p::Int) -> last(divide(num,den)(p))