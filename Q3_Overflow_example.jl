using Pkg
Pkg.activate(".")

include("poly_factorization_project.jl")


psparse=PolynomialSparse([Term(10,10)])#,Term(10,0)])
pBI=PolynomialSparseBI([BTerm(BigInt(10),10)])#,BTerm(BigInt(10),0)])

println("PolynomialSparse")
println("(",psparse,")","^30")
println(psparse^30)

println("\nPolynomialSparseBI")
println("(",pBI,")","^30")
println(pBI^30)

#Clearly (10⋅x^10)^10 should always be a multiple of 10
#PolynomialSparse wraps around giving an incorrect expression

#This becomes apparent for modulo operations 
println("\nPolynomialSparse")
println(mod(pBI^30,10))
println("\nPolynomialSparseBI")
println(mod(psparse^30,10))

#Hence this would interfere with more complicated function such as GCD that rely on mod and multiplication