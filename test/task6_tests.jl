

prime_samples=[2,3,5,7,11,13,23,29,31,37,41,43,47,53,59]

function refactored_pow_mod_test(;N::Int = 10^2, N_prods::Int = 5, seed::Int = 0)
    Random.seed!(seed)
    prime=rand(prime_samples)
    for _ in 1:N
        for i in 1:N_prods
        p1 = rand(PolynomialSparseBI)
        prod = pow_mod(p1,i,prime)
        @assert leading(prod) == mod(leading(p1)^i,prime) || mod(leading(p1)^i,prime) == BTerm(BigInt(0),0)
        end
    end
    println("refactored_pow_mod_test - PASSED")
end


function refactored_power_test(;N::Int = 10^2, N_prods::Int = 5, seed::Int = 0)
    Random.seed!(seed)
    prime=rand(prime_samples)
    for _ in 1:N
        for i in 1:N_prods

            p1 = rand(PolynomialSparse)
            p2 = PolynomialModP(p1,prime)
            prod=p2^i
            @assert leading(prod) == mod(leading(p1)^i,prime) || mod(leading(p1)^i,prime) == Term(0,0)
        end
    end
    println("refactored_power_test - PASSED")
end