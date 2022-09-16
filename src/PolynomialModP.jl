struct PolynomialModP
    poly::PolynomialSparse
    prime::Int
    function PolynomialModP(poly::PolynomialSparse,prime::Int)
        return(new(mod(poly,prime),prime))
    end
end

PolynomialModP(p::PolynomialModP, i::Int) = PolynomialModP(p,i)

+(pm::PolynomialModP,p::PolynomialModP) = PolynomialModP(mod((pm.poly+p.poly),pm.prime),pm.prime)

-(pm::PolynomialModP,p::PolynomialModP) = PolynomialModP(mod((pm.poly-p.poly),pm.prime),pm.prime)

*(pm::PolynomialModP,p::PolynomialModP) = PolynomialModP(mod((pm.poly*p.poly),pm.prime),pm.prime)

show(io::IO, p::PolynomialModP) = print(p.poly)

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