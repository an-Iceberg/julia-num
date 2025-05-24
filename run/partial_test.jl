using Format
using Printf

include("../src/partial.jl")

x, y = -1, 1
x⃗ = [x, y]

f(x, y) = 3x * y^3 + 10x^2 * y + 5y + 3y * sin(5x * y)
∂xf(x, y) = 3y^3 + 20x * y + 3y * cos(5x * y) * 5y
∂yf(x, y) = 3x * 3y^2 + 10x^2 + 5 + (3sin(5x * y) + 3y * cos(5x * y) * 5x)

function g(input)
  x, y = input
  return f(x, y)
end

@printf("%.2e\n", abs(∂xf(x, y) - ∂(g, x⃗, 1)))
@printf("%.2e\n", abs(∂yf(x, y) - ∂(g, x⃗, 2)))
