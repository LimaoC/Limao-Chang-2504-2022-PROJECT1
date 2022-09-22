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
integers_tests()

####
# Execute unit tests for polynomial operations
####
include("polynomial_product_test.jl")
polynomial_product_tests()

include("polynomial_derivative_test.jl")
polynomial_derivative_tests()

include("polynomial_ext_euclid_test.jl")
polynomial_ext_euclid_tests()

include("polynomial_division_test.jl")
polynomial_division_tests()

include("polynomial_overflow_test.jl")
polynomial_overflow_tests()

####
# Execute unit tests for polynomial factorization
####
include("factorization_test.jl")
factorization_tests()