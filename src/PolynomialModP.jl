struct PolynomialModP
    poly::PolynomialSparse
    prime::Int
    function PolynomialModP(poly::PolynomialSparse,prime::Int)
        return(new(mod(poly,prime),prime))
    end
end

prime(p::PolynomialModP) = p.prime

poly(p::PolynomialModP) = p.poly

PolynomialModP(p::PolynomialModP, i::Int) = PolynomialModP(poly(p),i)

+(pm::PolynomialModP,p::PolynomialModP) = PolynomialModP((pm.poly+p.poly),pm.prime)

-(pm::PolynomialModP,p::PolynomialModP) = PolynomialModP((pm.poly-p.poly),pm.prime)

*(pm::PolynomialModP,p::PolynomialModP) = PolynomialModP((pm.poly*p.poly),pm.prime)

*(n::Int,p::PolynomialModP) = PolynomialModP(n*p.poly,p.prime)

*(pm::PolynomialModP,p::PolynomialSparse) = PolynomialModP((pm.poly*p),pm.prime)

show(io::IO, p::PolynomialModP) = print(poly(p), " (mod ",prime(p),")")

degree(p::PolynomialModP) = degree(p.poly)


leading(p::PolynomialModP) = leading(p.poly)

iszero(p::PolynomialModP) = iszero(p.poly)

function divide(num::PolynomialModP, den::PolynomialModP)
        f, g = num, den
        p=num.prime
        degree(f) < degree(num) && return nothing 
        iszero(g) && throw(DivideError())
        q = PolynomialModP(PolynomialSparse(),p)
        prev_degree = degree(f)
        while degree(f) ≥ degree(g) 
            h = PolynomialModP(PolynomialSparse( (leading(f) ÷ leading(g))(p) ),p)  #syzergy 
            f = f - (h*g)
            q = q + h
            prev_degree == degree(f) && break
            prev_degree = degree(f)
        end
        @assert iszero( (num  - (q*g + f)))
        return q, f
end

÷(num::PolynomialModP, den::PolynomialModP)  = first(divide(num,den))

rem(num::PolynomialModP, den::PolynomialModP)  = last(divide(num,den))



function ^(p::PolynomialModP, n::Int)
    n < 0 && error("No negative power")
        #Creates a binary representation of n
        m=digits(n, base=2)
        ans=1
        w=p
        for i in m
            if i==1
                ans*=w
            end
            w*=w
        end
        return ans
    end


prim_part(p::PolynomialModP) =PolynomialModP((poly(p) ÷ content(poly(p)))(prime(p)),prime(p))




function extended_euclid_alg(a::PolynomialModP, b::PolynomialModP, prime::Int)
    old_r, r = PolynomialModP(poly(a),prime), PolynomialModP(poly(b),prime)
    old_s, s = PolynomialModP(one(PolynomialSparse),prime), PolynomialModP(zero(PolynomialSparse),prime)
    old_t, t = PolynomialModP(zero(PolynomialSparse),prime), PolynomialModP(one(PolynomialSparse),prime)

    while !iszero(PolynomialModP(poly(r),prime))
        q = first(divide(old_r, r))
        old_r, r = r, PolynomialModP(old_r - q*r, prime)
        old_s, s = s, PolynomialModP(old_s - q*s, prime)
        old_t, t = t, PolynomialModP(old_t - q*t, prime)
    end
    g, s, t = old_r, old_s, old_t

    @assert poly(PolynomialModP(s*a + t*b - g, prime)) == 0
    return g, s, t  
end

gcd(a::PolynomialModP, b::PolynomialModP, prime::Int) = extended_euclid_alg(a,b,prime) |> first



"""
Factoization methods
"""


function dd_split(f::PolynomialModP, d::Int,)::Vector{PolynomialModP}
    degree(f) == d && return [f]
    degree(f) == 0 && return []
    w = rand(PolynomialSparse, degree = d, monic = true)
    w=PolynomialModP(w,prime(f))
    n_power = (prime(f)^d-1) ÷ 2
    h=w^n_power - PolynomialModP(one(PolynomialSparse),prime(f))
    g = gcd(h,f,prime(f))
    ḡ = (f ÷ g) # g\bar + [TAB]
    return vcat(dd_split(g, d), dd_split(ḡ, d) )
end

function dd_factor(f::PolynomialModP)::Array{PolynomialModP}
    x = PolynomialModP(PolynomialSparse(Term(1,1)),prime(f))
    w = deepcopy(x)
    g = Array{PolynomialModP}(undef,degree(f)) #Array of polynomials indexed by degree

    #Looping over degrees
    for k in 1:degree(f)
        w = rem(w^prime(f), f)
        g[k] = gcd(w - x, f, prime(f)) 
        f = (f ÷ g[k])
    end


    #edge case for final factor
    f != PolynomialModP(one(PolynomialSparse),prime(f)) && push!(g,f)
    
    return g
end

function multiplicity(f::PolynomialModP, g::PolynomialModP, prime::Int)::Int
    degree(gcd(f, g, prime)) == 0 && return 0
    return 1 + multiplicity((f ÷ g), g, prime)
end

function expand_factorization(factorization::Vector{Tuple{PolynomialModP,Int}})::PolynomialModP
    length(factorization) == 1 && return first(factorization[1])^last(factorization[1])
    return *([first(tt)^last(tt) for tt in factorization]...)
end

function expand_factorization(factorization::Vector{Tuple{PolynomialSparse,Int}})::PolynomialSparse
    length(factorization) == 1 && return first(factorization[1])^last(factorization[1])
    return *([first(tt)^last(tt) for tt in factorization]...)
end

function factor(f::PolynomialModP)::Vector{Tuple{PolynomialModP,Int}}
    #Cantor Zassenhaus factorization

    degree(f) ≤ 1 && return [(f,1)]

    # make f primitive
    ff = prim_part(f)   
    # @show "after prim:", ff

     # make f square-free
    squares_poly = gcd(f, PolynomialModP(derivative(poly(ff)), prime(f)), prime(f)) 
    ff = (ff ÷ squares_poly)
    # @show "after square free:", ff

    # make f monic
    old_coeff = leading(ff).coeff
    ff = (poly(ff) ÷ old_coeff)(prime(f))      
    ff=PolynomialModP(ff,prime(f))
    # @show "after monic:", ff

    dds = dd_factor(ff)

    ret_val=Tuple{PolynomialModP,Int}[]

    for (k,dd) in enumerate(dds)
        sp = dd_split(dd, k)
        sp = map((p)->PolynomialModP((poly(p) ÷ leading(poly(p)).coeff)(prime(f)),prime(f)),sp) #makes the polynomials inside the list sp, monic
        for mp in sp
            push!(ret_val, (mp,multiplicity(f,mp,prime(f))))
        end
    end

    return ret_val
end


"""
tests
"""



"""
Test product of polynomials.
"""
function prod_test_polymodp(;N::Int = 10^3, N_prods::Int = 20, seed::Int = 0)
    Random.seed!(seed)
    primes=[3,5,7,11,13,17,23,29,31]
    for _ in 1:N
        prime1=rand(primes)
        p1 = PolynomialModP(rand(PolynomialSparse),prime1)
        p2 = PolynomialModP(rand(PolynomialSparse),prime1)
        prod = p1*p2
        @assert leading(prod) == mod(leading(p1)*leading(p2),prime1)
    end

    for _ in 1:N
        prime2=rand(primes)
        p_base = PolynomialModP(PolynomialSparse(Term(1,0)),prime2)
        for _ in 1:N_prods
            p = PolynomialModP(rand(PolynomialSparse),prime2)
            prod = p_base*p
            @assert leading(prod) == mod(leading(p_base)*leading(p),prime2)
            p_base = prod
        end
    end
    println("prod_test_polymodp - PASSED")
end


"""
Test the extended euclid algorithm for polynomials modulo p.
"""
function ext_euclid_test_poly(;prime::Int=101, N::Int = 10^3, seed::Int = 0)
    Random.seed!(seed)
    for _ in 1:N
        p1 = PolynomialModP(rand(PolynomialSparse),prime)
        p2 = PolynomialModP(rand(PolynomialSparse),prime)
        g, s, t = extended_euclid_alg(p1, p2, prime)
        @assert poly(s*p1 + t*p2 - g) == 0
    end
    println("ext_euclid_test_polymodp - PASSED")
end

p=PolynomialModP(rand(PolynomialSparse),101)
println(p)
@show factor(p)