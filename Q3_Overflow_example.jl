using Pkg
Pkg.activate(".")

include("poly_factorization_project.jl")


psparse=PolynomialSparse([Term(10,10)])#,Term(10,0)])
pBI=PolynomialSparseBI([BTerm(BigInt(10),10)])#,BTerm(BigInt(10),0)])

println("PolynomialSparse")
println("(",psparse,")","^30 =")
println(psparse^30)

println("\nPolynomialSparseBI")
println("(",pBI,")","^30 =")
println(pBI^30)

#Clearly (10⋅x^10)^10 should always be a multiple of 10

psparse2=PolynomialSparse([Term(10,10),Term(10,5),Term(10,0)])#,Term(10,0)])
pBI2=PolynomialSparseBI([BTerm(BigInt(10),10),BTerm(BigInt(10),5),BTerm(BigInt(10),0)])#,BTerm(BigInt(10),0)])

println("PolynomialSparse")
println("(",psparse2,")","^30 =")
println(psparse2^30)

println("\nPolynomialSparseBI")
println("(",pBI2,")","^30 =")
println(pBI2^30)

#And (10⋅x^10 + 10⋅x^5 + 10)^10 should always have coefficients with multiples of ten 

#PolynomialSparse wraps around giving an incorrect expression

#This becomes apparent for modulo operations 
println("\nPolynomialSparse")
println(mod(pBI^30,10))
println("\nPolynomialSparseBI")
println(mod(psparse^30,10))

#Hence this would interfere with more complicated functions such as GCD that rely on mod and multiplication