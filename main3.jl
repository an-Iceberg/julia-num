using ForwardDiff, Printf

D = ForwardDiff.jacobian
∇ = ForwardDiff.gradient
d = ForwardDiff.derivative

d2(f, x) = d(x -> d(f, x), x)
d3(f, x) = d2(x -> d(f, x), x)

f(x) = sin(x) - cos(x)
d1_f(x) = cos(x) + sin(x)
d2_f(x) = -sin(x) + cos(x)
d3_f(x) = -cos(x) - sin(x)

@printf("df = %.4f\n", d(f, 2))
@printf("d2f = %.4f\n", d2(f, 2))
@printf("d3f = %.4f\n", d3(f, 2))
println()
@printf("df ε = %.2e\n", abs(d1_f(2) - d(f, 2)))
@printf("d2f ε = %.2e\n", abs(d2_f(2) - d2(f, 2)))
@printf("d3f ε = %.2e\n", abs(d3_f(2) - d3(f, 2)))
println()

x = [-1, 1]

function f(input)
  x, y = input
  return 3x * y^3 + 10x^2 * y + 5y + 3y * sin(5x * y)
end

println()
display(∇(f, x))

function f(input)
  x, y = input
  return [
    3x * y^3 + 10x^2 * y + 5y + 3y * sin(5x * y),
    0,
  ]
end

println()
display(D(f, x))
