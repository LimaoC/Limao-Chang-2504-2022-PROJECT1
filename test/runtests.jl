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

include("integers_test.jl")
include("polynomial_product_test.jl")
include("polynomial_derivative_test.jl")
include("polynomial_ext_euclid_test.jl")
include("polynomial_division_test.jl")
include("polynomial_sparse_test.jl")
include("polynomial_overflow_test.jl")
include("factorization_test.jl")


# Integer unit tests
integers_tests()
# Polynomial operation unit tests
polynomial_product_tests()
polynomial_derivative_tests()
polynomial_ext_euclid_tests()
polynomial_division_tests()
# PolynomialSparse specific tests
polynomial_sparse_tests()
# PolynomialSparseBI specific tests
polynomial_overflow_tests()
# Polynomial factorization tests
factorization_tests()