using Format, ForwardDiff

D = ForwardDiff.jacobian
∇ = ForwardDiff.gradient
d = ForwardDiff.derivative

f(x) = sin(x) - cos(x)

# display(ForwardDiff.jacobian(f, [2]))
# display(D(f, [2]))
# println()
# display(ForwardDiff.gradient(f, 2))
# display(∇(f, 2))
# println()
printfmtln("{:.4f}", ForwardDiff.derivative(f, 2))
printfmtln("{:.4f}", d(f, 2))
printfmtln("{:.4f}", d(x -> d(f, x), 2))
printfmtln("{:.4f}", d(x -> d(x -> d(f, x), x), 2))

x = [-1, 1]

function f(input)
  x, y = input
  return 3x * y^3 + 10x^2 * y + 5y + 3y * sin(5x * y)
end

println()
display(∇(f, x))

function f(input)
  x, y = input
  return [3x * y^3 + 10x^2 * y + 5y + 3y * sin(5x * y), 0]
end

println()
display(D(f, x))
