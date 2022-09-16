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
        if !iszero(p.terms[t.degree + 1]) #+1 is due to indexing
            p.terms[t.degree + 1] += t
        else
            p.terms[t.degree + 1] = t
        end
    end

    return trim!(p)
end

function +(p::PolynomialSparse, t::Term)
    p = deepcopy(p)
    if degree(t) ∉ degree.(p.terms)
        push!(p, t)
    else
            loc=1
            for i in p
                if degree(i)==degree(t)
                    p.terms[loc]=Term(coeff(i)+coeff(t),degree(t))
                    break
                else loc+=1
                end
            end
    end

    return PolynomialSparse(p.terms)
end

+(t::Term, p::Union{PolynomialDense, PolynomialSparse}) = p + t

function +(p::PolynomialSparseBI, t::BTerm)
    p = deepcopy(p)
    if degree(t) ∉ degree.(p.terms)
        push!(p, t)
    else
            loc=1
            for i in p
                if degree(i)==degree(t)
                    p.terms[loc]=BTerm(coeff(i)+coeff(t),degree(t))
                    break
                else loc+=1
                end
            end
    end

    return PolynomialSparseBI(p.terms)
end

+(t::Term, p::Union{PolynomialDense, PolynomialSparse, PolynomialSparseBI}) = p + t



"""
Add two polynomials.
"""
function +(p1::PolynomialDense, p2::PolynomialDense)::PolynomialDense
    p = deepcopy(p1)
    for t in p2
        p += t
    end
    return p
end

function +(p1::PolynomialSparse, p2::PolynomialSparse)::PolynomialSparse
    p = deepcopy(p1)
    for t in p2

        p += t
    end
    return p
end

function +(p1::PolynomialSparseBI, p2::PolynomialSparseBI)::PolynomialSparseBI
    p = deepcopy(p1)
    for t in p2

        p += t
    end
    return p
end

"""
Add a polynomial and an integer.
"""
+(p::Union{PolynomialDense, PolynomialSparse}, n::Int) = p + Term(n,0)
+(n::Int, p::Union{PolynomialDense, PolynomialSparse}) = p + Term(n,0)

+(p::PolynomialSparseBI, n::Int) = p + Term(BigInt(n),0)
+(n::Int, p::PolynomialSparseBI) = p + Term(BigInt(n),0)