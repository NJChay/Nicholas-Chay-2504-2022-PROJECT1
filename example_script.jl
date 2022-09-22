using Pkg
Pkg.activate(".")

include("poly_factorization_project.jl")

x = x_poly()
p1 = 2x^3 + 4x^2 - 3x
p2 = 2x^4 - 4x^2 - 3x + 3

@show p1*p2
@show p1^3

@show derivative(p1*p2)
@show derivative(p1)*p2 + p1*derivative(p2);

prime1 = 17
p = mod((7x^3 + 2x^2 + 8x + 1)*(x^2+x+1),prime1)
println("Will factor this polynomial (mod $prime1): ", p)
factorization = factor(p,prime1)
println("Here is the factorization: ", factorization)

pr = mod(expand_factorization(factorization),prime1)
println("Reconstructing: ", pr)

x = x_poly()
p1 = 2x^2+3x +(-5) #Note we need +(-5)... try later without... how to fix?
p2 = -3x^2 - 4x +6
p1+p2

extended_euclid_alg(p1*p2,p2,101)