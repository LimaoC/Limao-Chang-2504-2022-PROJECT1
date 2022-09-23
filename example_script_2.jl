using Pkg
Pkg.activate(".")

include("poly_factorization_project.jl")

x = x_poly(PolynomialDense)
dense_polynomial = x^5 + x
prime = 17
p1 = x^2 + 5
p2 = x^3 + 12x + 6
x = x_poly(PolynomialSparse)
sparse_polynomial = x^5 + x
factorization = factor(p1 * p2, prime)
pr = mod(expand_factorization(factorization), prime)

println("We can construct polynomials using the basic `x` polynomial, like so:")
println("First we define `x = x_poly(PolynomialDense)`. Then we can construct " *
        "polynomials using intuitive addition, subtraction, and multiplication " *
        "operations like so:")

println(2x)
println(3x^2)
println(5x^4 - 2x + 6)

println("\nWe can also perform polynomial on polynomial arithmetic. Given two " *
        "polynomials p1 and p2, we can perform:")
println("addition:")
@show p1 + p2
println("subtraction:")
@show p1 - p2
println("multiplication:")
@show p1 * p2
println("division modulo p (say p = 17):")
println("p2 ÷ p1 = $((p2 ÷ p1)(prime))")
println("differentiation:")
@show derivative(p1 * p2)
@show derivative(p1) * p2 + p1 * derivative(p2);

println("\nFor polynomials of large degree, we may consider using a sparse " *
        "polynomial representation over a dense polynomial representation. Take the " *
        "polynomial x^5 + x for example:")
println("Dense polynomial terms: $(dense_polynomial.terms)")
println("Sparse polynomial terms: $(sparse_polynomial.terms)")
println("Sparse polynomials only keep non-zero terms to be more space efficient - we " *
        "don't need to store any unnecessary zeroes.")

println("\nAnother polynomial operation we can perform is factoring. Taking the same " *
        "p1 and p2, we can factor the polynomial modulo 17 like so:")
println(factorization)

println("And now building p1 * p2 back up again using these factors:")
println(pr)
