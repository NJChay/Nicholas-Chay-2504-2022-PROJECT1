function CRT(v::Vector{Any})
    z=0
    total=1
    for i in v
        total=total*i[2]
    end
    for i in v
        b=total ÷ i[2]
        z+= b*i[1]*int_inverse_mod(b,i[2])
    end
    while z > total
        z=z-total
    end
    return smod.(z,total)
end

smod(a::Int, m::Int) = mod(a,m) <= m//2 ? mod(a,m) : mod(a,m) - m

function CRTm(a::PolynomialSparseBI,b::PolynomialSparseBI)
#random primes
    primes=[599,601,409]

    modpolys=[]
    sorted=[]
    #Convert components to a list of terms mod each prime
    for i in primes
        newp=PolynomialModP(PolynomialSparse(a),i)*PolynomialModP(PolynomialSparse(b),i)
        push!(modpolys,(poly(newp).terms))
    end

    #Fills in 0 for unspecified terms less than the maximum degree

    size=maximum(degree.(modpolys[1]))+1
    for p in modpolys 
        terms = [zero(Term) for i in 0:size-1]

        for (n,t) in enumerate(p)
            terms[t.degree+1] = t 
        end
        push!(sorted,reverse(terms))

    end

    #Resorts coefficients and groups them with a prime 

    sorted=[coeff.(i) for i in sorted]
    resorted=[]
    
    for i in 1:size
        new_row=[]
        for (ii,t) in enumerate(sorted)
            push!(new_row,[t[i],primes[ii]])
        end
        push!(resorted,new_row)
    end

    #Constructs new polynomial from the CRT result
    Bterms::Vector{BTerm}=[]
    i=0

    for t in reverse(CRT.(resorted))
        push!(Bterms,BTerm(BigInt(t),i))
        i+=1
    end

    return PolynomialSparseBI(Bterms)

end

function CRT_multiplication_test(N::Int = 10^3, seed::Int = 0)
    Random.seed!(seed)
    for _ in 1:N
        p1 = rand(PolynomialSparseBI)
        p2 = rand(PolynomialSparseBI)
        prod = CRTm(p1,p2)
        @assert leading(prod) == leading(p1)*leading(p2)
    end
    println("CRT_multiplication_test - PASSED")
end