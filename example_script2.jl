using Pkg
Pkg.activate(".")

include("poly_factorization_project.jl")

println("Example Polynomials")
p1=(rand(PolynomialSparse))
p2=(rand(PolynomialSparse))
p3=(rand(PolynomialSparse))
p4=(rand(PolynomialSparse))

println("p1 = ", p1)
println("p2 = ", p2)
println("p3 = ", p3)
println("p4 = ", p4)


println("\nModulo")
println("p1 modulo 2: p1=", mod(p1,2))
println("p2 modulo 3: p2=", mod(p2,3))
println("p3 modulo 5: p3=", mod(p3,5))
println("p4 modulo 7: p4=", mod(p4,7))

p1=mod(p1,2)
p2=mod(p2,3)
p3=mod(p3,5)
p4=mod(p4,7)

println("\nAddition")
println("p1+p2 =",p1+p2)
println("p3+p4 =",p3+p4)

println("\nGCD")
 println("GCD of p1,p2 and 6 =",gcd(p1,p2,6))
 println("GCD of p3,p4 and 20 =",gcd(p3,p4,20))

println("\nMultiplication")
 println("p1*p2 =", p1*p2)
 println("p3*p4 =", p3*p4)

 println("\nExponents")
println("p1^4 =", p1^4)
println("p1^8 =", p2^8)

println("\nDerivative")

println("\nFactorization")
 println("factors of p1 with prime 3:",factor(p1,3))
 println("fatcors of p3 with prime 5:", factor(p3,5))

println("\nReconstruction")

