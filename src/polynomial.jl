#############################################################################
#############################################################################
#
# This file defines the polynomial type with several operations 
#                                                                               
#############################################################################
#############################################################################

####################################
# Polynomial type and construction #
####################################


"""
A Polynomial type - designed to be for polynomials with integer coefficients.
"""

struct Polynomial
end



struct PolynomialDense

    #A zero packed vector of terms
    #Terms are assumed to be in order with first term having degree 0, second degree 1, and so fourth
    #until the degree of the polynomial. The leading term (i.e. last) is assumed to be non-zero except 
    #for the zero polynomial where the vector is of length 1.
    #Note: at positions where the coefficient is 0, the power of the term is also 0 (this is how the Term type is designed)
    terms::Vector{Term}   
    
    #Inner constructor of 0 polynomial
    PolynomialDense() = new([zero(Term)])

    #Inner constructor of polynomial based on arbitrary list of terms
    function PolynomialDense(vt::Vector{Term})

        #Filter the vector so that there is not more than a single zero term
        vt = filter((t)->!iszero(t), vt)
        if isempty(vt)
            vt = [zero(Term)]
        end

        max_degree = maximum((t)->t.degree, vt)
        terms = [zero(Term) for i in 0:max_degree] #First set all terms with zeros

        #now update based on the input terms
        for t in vt
            terms[t.degree + 1] = t #+1 accounts for 1-indexing
        end
        return new(terms)
    end
end



struct PolynomialSparse
    terms::Vector{Term}  
    
    #Inner constructor of 0 polynomial
    PolynomialSparse() = new([zero(Term)])

    #Inner constructor of polynomial based on arbitrary list of terms
    function PolynomialSparse(vt::Vector{Term})

        #Filter the vector so that there is not more than a single zero term
        vt = filter((t)->!iszero(t), vt)
        if isempty(vt)
            vt = [zero(Term)]
        end
        vt = sort(vt, by=degree, rev=true)

        return new(vt)
    end
end



struct PolynomialSparseBI
    terms::Vector{BTerm}  
    
    #Inner constructor of 0 polynomial
    PolynomialSparseBI() = new([zero(BTerm)])

    #Inner constructor of polynomial based on arbitrary list of terms
    function PolynomialSparseBI(vt::Vector{BTerm})

        #Filter the vector so that there is not more than a single zero term
        vt = filter((t)->!iszero(t), vt)
        if isempty(vt)
            vt = [zero(BTerm)]
        end
        vt = sort(vt, by=degree, rev=true)

        return new(vt)
    end
end



"""
This function maintains the invariant of the Polynomial type so that there are no zero terms beyond the highest
non-zero term.
"""
function trim!(p::Union{PolynomialDense,PolynomialSparse, PolynomialSparseBI})
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
Construct a polynomial with a single term.
"""
PolynomialDense(t::Term) = PolynomialDense([t])
PolynomialSparse(t::Term) = PolynomialSparse([t])
PolynomialSparseBI(t::BTerm) = PolynomialSparseBI([t])


"""
Construct a polynomial of the form x^p-x.
"""
cyclotonic_polynomial(p::Int) = PolynomialDense([Term(1,p), Term(-1,0)])


"""
Construct a polynomial of the form x-n.
"""
linear_monic_polynomial(n::Int) = PolynomialDense([Term(1,1), Term(-n,0)])

"""
Construct a polynomial of the form x.
"""
x_poly() = PolynomialDense(Term(1,1))

"""
Creates the zero polynomial.
"""
zero(::Type{PolynomialDense})::PolynomialDense = PolynomialDense()
zero(::Type{PolynomialSparse})::PolynomialSparse = PolynomialSparse()
zero(::Type{PolynomialSparseBI})::PolynomialSparseBI = PolynomialSparseBI()

"""
Creates the unit polynomial.
"""
one(::Type{PolynomialDense})::PolynomialDense = PolynomialDense(one(Term))
one(p::PolynomialDense) = one(typeof(p))

one(::Type{PolynomialSparse})::PolynomialSparse = PolynomialSparse(one(Term))
one(p::PolynomialSparse) = one(typeof(p))

one(::Type{PolynomialSparseBI})::PolynomialSparseBI = PolynomialSparseBI(one(BTerm))
one(p::PolynomialSparseBI) = one(typeof(p))
"""
Generates a random polynomial.
"""
function rand(::Type{PolynomialDense} ; 
                degree::Int = -1, 
                terms::Int = -1, 
                max_coeff::Int = 100, 
                mean_degree::Float64 = 5.0,
                prob_term::Float64  = 0.7,
                monic = false,
                condition = (p)->true)
        
    while true 
        _degree = degree == -1 ? rand(Poisson(mean_degree)) : degree
        _terms = terms == -1 ? rand(Binomial(_degree,prob_term)) : terms
        degrees = vcat(sort(sample(0:_degree-1,_terms,replace = false)),_degree)
        coeffs = rand(1:max_coeff,_terms+1)
        monic && (coeffs[end] = 1)
        p = PolynomialDense( [Term(coeffs[i],degrees[i]) for i in 1:length(degrees)] )
        condition(p) && return p
    end
end

function rand(::Type{PolynomialSparse} ; 
    degree::Int = -1, 
    terms::Int = -1, 
    max_coeff::Int = 100, 
    mean_degree::Float64 = 5.0,
    prob_term::Float64  = 0.7,
    monic = false,
    condition = (p)->true)
    p=nothing
    while true 
    _degree = degree == -1 ? rand(Poisson(mean_degree)) : degree
    _terms = terms == -1 ? rand(Binomial(_degree,prob_term)) : terms
    degrees = vcat(sort(sample(0:_degree-1,_terms,replace = false)),_degree)
    coeffs = rand(1:max_coeff,_terms+1)
    monic && (coeffs[end] = 1)
    p = PolynomialSparse( [Term(coeffs[i],degrees[i]) for i in 1:length(degrees)] )
    condition(p) && break
    end
    return p
end

function rand(::Type{PolynomialSparseBI} ; 
    degree::Int = -1, 
    terms::Int = -1, 
    max_coeff::Int = 100, 
    mean_degree::Float64 = 5.0,
    prob_term::Float64  = 0.7,
    monic = false,
    condition = (p)->true)
    p=nothing
    while true 
    _degree = degree == -1 ? rand(Poisson(mean_degree)) : degree
    _terms = terms == -1 ? rand(Binomial(_degree,prob_term)) : terms
    degrees = vcat(sort(sample(0:_degree-1,_terms,replace = false)),_degree)
    coeffs = rand(1:max_coeff,_terms+1)
    monic && (coeffs[end] = 1)
    p = PolynomialSparseBI( [BTerm(BigInt(coeffs[i]),degrees[i]) for i in 1:length(degrees)] )
    condition(p) && break
    end
    return p
end


###########
# Display #
###########

"""
Show a polynomial.
"""
function show(io::IO, p::Union{PolynomialDense, PolynomialSparse}) 
    if iszero(p)
        print(io,"0")
    else
        if typeof(p)==PolynomialSparse
            true_terms=[term for term in p.terms if term.coeff!=0]
        else
            true_terms=reverse([term for term in p.terms if term.coeff!=0])
        end
        n = length(true_terms)
        for (i,t) in enumerate(true_terms)
            if !iszero(t) && n!=1
                abs_term=Term(abs(t.coeff),t.degree)
                if i!=1 && i<n && true_terms[i+1].coeff>=0
                    print(io, abs_term, i != n ? " + " : "")
                    
                elseif  i==1 && true_terms[i+1].coeff<0
                    print( t, i != n ? " - " : "")

                elseif i==1
                    print( t, i != n ? " + " : "")
                else
                    print(abs_term, i != n ? " - " : "")
                end
            elseif !iszero(t)
                print(t)
            end
        end
    end
end

function show(io::IO, p::PolynomialSparseBI) 
    if iszero(p)
        print(io,"0")
    else
        true_terms=[term for term in p.terms if term.coeff!=0]
        n = length(true_terms)
        for (i,t) in enumerate(true_terms)
            if !iszero(t) && n!=1
                abs_term=BTerm(abs(t.coeff),t.degree)
                if i!=1 && i<n && true_terms[i+1].coeff>=0
                    print(io, abs_term, i != n ? " + " : "")
                    
                elseif  i==1 && true_terms[i+1].coeff<0
                    print( t, i != n ? " - " : "")

                elseif i==1
                    print( t, i != n ? " + " : "")
                else
                    print(abs_term, i != n ? " - " : "")
                end
            elseif !iszero(t)
                print(t)
            end
        end
    end
end
##############################################
# Iteration over the terms of the polynomial #
##############################################

"""
Allows to do iteration over the non-zero terms of the polynomial. This implements the iteration interface.
"""
iterate(p::Union{PolynomialDense, PolynomialSparse, PolynomialSparseBI}, state=1) = iterate(p.terms, state)


##############################
# Queries about a polynomial #
##############################

"""
The number of terms of the polynomial.
"""
length(p::Union{PolynomialDense, PolynomialSparse, PolynomialSparseBI}) = length(p.terms) 
"""
The leading term of the polynomial.
"""
leading(p::PolynomialDense)::Term = isempty(p.terms) ? zero(Term) : last(p.terms)  
leading(p::PolynomialSparse)::Term = isempty(p.terms) ? zero(Term) : (p.terms)[1]
leading(p::PolynomialSparseBI)::BTerm = isempty(p.terms) ? zero(BTerm) : (p.terms)[1]
"""
Returns the coefficients of the polynomial.
"""
coeffs(p::Union{PolynomialDense, PolynomialSparse})::Vector{Int} = [t.coeff for t in p]
coeffs(p::Union{PolynomialSparseBI})::Vector{BigInt} = [t.coeff for t in p]
"""
The degree of the polynomial.
"""
degree(p::Union{PolynomialDense, PolynomialSparse, PolynomialSparseBI})::Int = leading(p).degree 

"""
The content of the polynomial is the GCD of its coefficients.
"""
content(p::Union{PolynomialDense, PolynomialSparse, PolynomialSparseBI})::Int = euclid_alg(coeffs(p))

"""
Evaluate the polynomial at a point `x`.
"""
evaluate(f::Union{PolynomialDense, PolynomialSparse, PolynomialSparseBI}, x::T) where T <: Number = sum(evaluate(t,x) for t in f)

################################
# Pushing and popping of terms #
################################

"""
Push a new term into the polynomial.
"""
#Note that ideally this would throw and error if pushing another term of degree that is already in the polynomial
function push!(p::PolynomialDense, t::Term) 
    if t.degree <= degree(p)
        p.terms[t.degree + 1] = t
    else
        append!(p.terms, zeros(Term, t.degree - degree(p)-1))
        push!(p.terms, t)
    end
    return p        
end

function push!(p::PolynomialSparse, t::Term) 
    p=PolynomialSparse(push!(p.terms,t))
    return p       
end

function push!(p::PolynomialSparseBI, t::BTerm) 
    p=PolynomialSparseBI(push!(p.terms,t))
    return p       
end


"""
Pop the leading term out of the polynomial. When polynomial is 0, keep popping out 0.
"""
function pop!(p::Union{PolynomialDense, PolynomialSparse})::Term 
    popped_term = pop!(p.terms) #last element popped is leading coefficient

    while !isempty(p.terms) && iszero(last(p.terms))
        pop!(p.terms)
    end

    if isempty(p.terms)
        push!(p.terms, zero(Term))
    end

    return popped_term
end

function pop!(p::PolynomialSparseBI)::BTerm 
    popped_term = pop!(p.terms) #last element popped is leading coefficient

    while !isempty(p.terms) && iszero(last(p.terms))
        pop!(p.terms)
    end

    if isempty(p.terms)
        push!(p.terms, zero(BTerm))
    end

    return popped_term
end

"""
Check if the polynomial is zero.
"""
iszero(p::Union{PolynomialDense,PolynomialSparse})::Bool = p.terms == [Term(0,0)]
iszero(p::PolynomialSparseBI)::Bool = length(p.terms)==1 && p.terms[1].coeff==0

#################################################################
# Transformation of the polynomial to create another polynomial #
#################################################################

"""
The negative of a polynomial.
"""
-(p::PolynomialDense) = PolynomialDense(map((pt)->-pt, p.terms))
-(p::PolynomialSparse) = PolynomialSparse(map((pt)->-pt, p.terms))
-(p::PolynomialSparseBI) = PolynomialSparseBI(map((pt)->-pt, p.terms))

"""
Create a new polynomial which is the derivative of the polynomial.
"""
function derivative(p::Union{PolynomialDense,PolynomialSparse})
    if typeof(p)==PolynomialSparse
        der_p = PolynomialSparse()
        Dense=false
    else
        der_p=PolynomialDense()
        Dense=true
    end
    for term in p
        der_term = derivative(term)
        !iszero(der_term) && push!(der_p,der_term)
    end
    if Dense==true
        return trim!(der_p)
    else
        return PolynomialSparse(trim!(der_p).terms)
    end
end

function derivative(p::PolynomialSparseBI)
    der_p = PolynomialSparseBI()
        
    for term in p
        der_term = derivative(term)
        !iszero(der_term) && push!(der_p,der_term)
    end
        return PolynomialSparseBI(trim!(der_p).terms)
end

"""
The prim part (multiply a polynomial by the inverse of its content).
"""
prim_part(p::Union{PolynomialDense, PolynomialSparse, PolynomialSparseBI}) = p ÷ content(p)


"""
A square free polynomial.
"""
square_free(p::Union{PolynomialDense, PolynomialSparse, PolynomialSparseBI}, prime::Int) = (p ÷ gcd(p,derivative(p),prime))(prime)

#################################
# Queries about two polynomials #
#################################

"""
Check if two polynomials are the same
"""
==(p1::Union{PolynomialDense, PolynomialSparse, PolynomialSparseBI}, p2::Union{PolynomialDense, PolynomialSparse, PolynomialSparseBI})::Bool = p1.terms == p2.terms


"""
Check if a polynomial is equal to 0.
"""
#Note that in principle there is a problem here. E.g The polynomial 3 will return true to equalling the integer 2.
==(p::Union{PolynomialDense,PolynomialSparse, PolynomialSparseBI}, n::T) where T <: Real = iszero(p) == iszero(n)

##################################################################
# Operations with two objects where at least one is a polynomial #
##################################################################

"""
Subtraction of two polynomials.
"""
-(p1::PolynomialSparse, p2::PolynomialSparse) = p1 + (-p2)
-(p1::PolynomialDense, p2::PolynomialDense) = p1 + (-p2)
-(p1::PolynomialSparseBI, p2::PolynomialSparseBI) = p1 + (-p2)

"""
Multiplication of polynomial and term.
"""
*(t::Term, p1::PolynomialDense) = iszero(t) ? PolynomialDense() : PolynomialDense(map((pt)->t*pt, p1.terms))
*(p1::PolynomialDense, t::Term)::PolynomialDense = t*p1

*(t::Term, p1::PolynomialSparse) = iszero(t) ? PolynomialSparse() : PolynomialSparse(map((pt)->t*pt, p1.terms))
*(p1::PolynomialSparse, t::Term)::PolynomialSparse = t*p1

*(t::BTerm, p1::PolynomialSparseBI) = iszero(t) ? PolynomialSparseBI() : PolynomialSparseBI(map((pt)->t*pt, p1.terms))
*(p1::PolynomialSparseBI, t::BTerm)::PolynomialSparseBI = t*p1
"""
Multiplication of polynomial and an integer.
"""
*(n::Int, p::PolynomialDense)::PolynomialDense = p*Term(n,0)
*(p::PolynomialDense, n::Int)::PolynomialDense = n*p

*(n::Int, p::PolynomialSparse)::PolynomialSparse = p*Term(n,0)
*(p::PolynomialSparse, n::Int)::PolynomialSparse = n*p

*(n::Int, p::PolynomialSparseBI)::PolynomialSparseBI = p*BTerm(BigInt(n),0)
*(p::PolynomialSparseBI, n::Int)::PolynomialSparseBI = n*p

*(n::BigInt, p::PolynomialSparseBI)::PolynomialSparseBI = p*BTerm(n,0)
*(p::PolynomialSparseBI, n::BigInt)::PolynomialSparseBI = n*p
"""
Integer division of a polynomial by an integer.

Warning this may not make sense if n does not divide all the coefficients of p.
"""
÷(p::PolynomialDense, n::Int) = (prime)->PolynomialDense(map((pt)->((pt ÷ n)(prime)), p.terms))
÷(p::PolynomialSparse, n::Int) = (prime)->PolynomialSparse(map((pt)->((pt ÷ n)(prime)), p.terms))
÷(p::PolynomialSparseBI, n::Int) = (prime)->PolynomialSparseBI(map((pt)->((pt ÷ n)(prime)), p.terms))
÷(p::PolynomialSparseBI, n::BigInt) = (prime)->PolynomialSparseBI(map((pt)->((pt ÷ n)(prime)), p.terms))

"""
Take the mod of a polynomial with an integer.
"""
function mod(f::PolynomialDense, p::Int)
    f_out = deepcopy(f)
    for i in 1:length(f_out.terms)
        f_out.terms[i] = mod(f_out.terms[i], p)
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


function mod(f::PolynomialSparse, p::Int)
    f_out = deepcopy(f)
    for i in 1:length(f_out.terms)
        f_out.terms[i] = mod(f_out.terms[i], p)
    end
    return PolynomialSparse(f_out.terms)
        
    # p_out = Polynomial()
    # for t in f
    #     new_term = mod(t, p)
    #     @show new_term
    #     push!(p_out, new_term)
    # end
    # return p_out
end

function mod(f::PolynomialSparseBI, p::Int)
    f_out = deepcopy(f)
    for i in 1:length(f_out.terms)
        f_out.terms[i] = mod(f_out.terms[i], p)
    end
    return PolynomialSparseBI(f_out.terms)
        
    # p_out = Polynomial()
    # for t in f
    #     new_term = mod(t, p)
    #     @show new_term
    #     push!(p_out, new_term)
    # end
    # return p_out
end

"""
Power of a polynomial mod prime.
"""
function pow_mod(p::Union{PolynomialDense, PolynomialSparse, PolynomialSparseBI}, n::Int, prime::Int)
    n < 0 && error("No negative power")
    out = one(p)
    for _ in 1:n
        out *= p
        out = mod(out, prime)
    end
    return out
end