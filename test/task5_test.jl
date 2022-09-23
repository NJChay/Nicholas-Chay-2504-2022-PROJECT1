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