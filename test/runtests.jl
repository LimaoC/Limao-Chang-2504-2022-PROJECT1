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
prod_test_poly_dense()
prod_test_poly_sparse()
prod_test_poly_sparse_bi()
prod_derivative_test_poly_dense()
prod_derivative_test_poly_sparse()
prod_derivative_test_poly_sparse_bi()
ext_euclid_test_poly_dense()
ext_euclid_test_poly_sparse()
ext_euclid_test_poly_sparse_bi()
division_test_poly_dense()
division_test_poly_sparse()
division_test_poly_sparse_bi()

####
# Execute unit tests for polynomial factorization
####
include("factorization_test.jl")
factor_test_poly_dense()
factor_test_poly_sparse()
factor_test_poly_sparse_bi()