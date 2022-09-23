#############################################################################
#############################################################################
#
# A script that runs all unit tests in the project.
#                                                                               
#############################################################################
#############################################################################

using Pkg
Pkg.activate(".")

include("../poly_factorization_project.jl")

# ~5 minutes for all tests


####
# Execute unit tests for integers
###
include("integers_test.jl")
test_euclid_ints()
test_ext_euclid_ints()

####
# Execute unit tests for polynomials
####
include("polynomials_test.jl")
include("polynomials_test _sparse.jl")
include("polynomials_test _sparse_BI.jl")
prod_test_poly()
prod_test_poly_sparse()
prod_test_poly_BI()
prod_derivative_test_poly()
prod_derivative_test_poly_sparse()
prod_derivative_test_poly_BI()
ext_euclid_test_poly()
ext_euclid_test_poly_sparse()
ext_euclid_test_poly_BI()
division_test_poly()
division_test_poly_sparse()
division_test_poly_BI()

####
# Execute unit tests for polynomial factorization
####
include("factorization_test.jl")
include("factorization_test _sparse.jl")
include("factorization_test _sparse_BI.jl")
factor_test_poly_BI()
factor_test_poly_sparse()
factor_test_poly()

####
# Execute unit tests for CRT multiplication and refined power and pow_mod
####

include("task5_test.jl")
include("task6_tests.jl")
CRT_multiplication_test()
refactored_pow_mod_test()
refactored_power_test()


