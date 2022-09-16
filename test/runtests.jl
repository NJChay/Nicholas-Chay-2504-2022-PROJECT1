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

# ~4 minutes for all tests


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
factor_test_poly_BI()
factor_test_poly_sparse()
factor_test_poly()