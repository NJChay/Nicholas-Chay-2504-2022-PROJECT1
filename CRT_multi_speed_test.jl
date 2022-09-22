using BenchmarkTools


p1=rand(PolynomialSparseBI)

p1=p1^2

println("regular multiplication")
@btime p1*p1

println("\nCRT multiplication")
@btime CRTm(p1,p1)

p1=p1^2

println("\nregular multiplication)")
@btime p1*p1


println("\nCRT multiplication")
@btime CRTm(p1,p1)


p2=rand(PolynomialSparseBI)

p2=p2^2

println("\nregular multiplication")
@btime p2*p2

println("\nCRT multiplication")
@btime CRTm(p2,p2)

println("\n")