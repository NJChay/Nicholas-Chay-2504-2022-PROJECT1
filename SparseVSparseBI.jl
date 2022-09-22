using BenchmarkTools

sparse=rand(PolynomialSparse, degree=10)
sparseBI=rand(PolynomialSparseBI, degree=10)


println("sparse")
@btime sparse^2
println("\nsparseBI")
@btime sparseBI^2

println("\nsparse")
@btime sparse^3
println("\nsparseBI")
@btime sparseBI^3

println("\nsparse")
@btime sparse^4
println("\nsparseBI")
@btime sparseBI^4

println("\n")