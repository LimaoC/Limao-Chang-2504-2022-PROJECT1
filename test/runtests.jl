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
include("polynomial_mod_p_test.jl")
include("polynomial_power_test.jl")
include("factorization_test.jl")


println("--- Integer unit tests ---")
@time integers_tests()

# Polynomial operation unit tests
println("--- Polynomial operations unit tests ---")
@time polynomial_product_tests()
@time polynomial_derivative_tests()
@time polynomial_ext_euclid_tests()
@time polynomial_division_tests()

println("--- Sparse polynomial representation unit tests ---")
@time polynomial_sparse_tests()

println("--- BigInt polynomial overflow unit tests ---")
@time polynomial_overflow_tests()

println("--- PolynomialModP specific tests ---")
@time polynomial_mod_p_tests()

println("--- Polynomial power tests ---")
@time polynomial_power_tests()

println("--- Polynomial factorization tests ---")
@time factorization_tests()