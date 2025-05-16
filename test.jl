using Printf

module A

function d(f::Function, x, h::Real=0.001)
  # return ((1 / 12)f(x - 2h) - (2 / 3)f(x - h) + (2 / 3)f(x + h) - (1 / 12)f(x + 2h)) / h

  # Even this implementation seems to be more accurate than the direct formulas?
  return (-f(x - h) + f(x + h)) / 2h
end

function d2(f::Function, x, h::Real=0.001)
  return d(x -> d(f, x, h), x)
end

# Precision seems to deteriorate after here
function d3(f::Function, x, h::Real=0.001)
  return d(x -> d2(f, x, h), x)
end

function d4(f::Function, x, h::Real=0.001)
  return d(x -> d3(f, x, h), x)
end

function ∇(f::Function, x⃗, h::Real=1e-3)
  h⃗ = fill(h, length(x⃗))
  return (-0.5f(x⃗ - h⃗) + 0.5f(x⃗ + h⃗)) / h
end

function ∇2(f::Function, x⃗, h::Real=1e-3)
  return ∇(x -> ∇(f, x⃗, h), x⃗, h)
end

end

module B

function d4_4(f::Function, x, h::Real=1e-3)
  return (
    -(1 / 6)f(x - 3h) + 2f(x - 2h) - (13 / 2)f(x - h) + (28 / 3)f(x) - (13 / 2)f(x + h) +
    2f(x + 2h) - (1 / 6)f(x + 3h)
  ) / h^4
end

function d4_6(f::Function, x, h::Real=1e-3)
  return (
    (7 / 240)f(x - 4h) - 0.4f(x - 3h) + (169 / 60)f(x - 2h) - (122 / 15)f(x - h) + 11.375f(x) -
    (122 / 15)f(x + h) + (169 / 60)f(x + 2h) - 0.4f(x + 3h) + (7 / 240)f(x + 4h)
  ) / h^4
end

end

f1(x) = x^3
f2(x) = 1 / x
f3(x) = sqrt(x)
f4(x) = exp(x)
f5(x) = log(x)

d4_f1(x) = 0
d4_f2(x) = 24 / x^5
d4_f3(x) = -(15 / 16x^(7 / 2))
d4_f4(x) = exp(x)
d4_f5(x) = -(6 / x^4)

# Todo: scoring
A_score = 0
B_score = 0

for x in [1, 2, 5, 10, 20, 50, 100, 200, 500]
  println("  x = $x")
  @printf("A f1 ε = %.2e\n", abs(d4_f1(x) - A.d4(f1, x)))
  @printf("B f1 ε = %.2e\n", abs(d4_f1(x) - B.d4_4(f1, x)))
  @printf("B f1 ε = %.2e\n", abs(d4_f1(x) - B.d4_6(f1, x)))
  println()
  @printf("A f2 ε = %.2e\n", abs(d4_f2(x) - A.d4(f2, x)))
  @printf("B f2 ε = %.2e\n", abs(d4_f2(x) - B.d4_4(f2, x)))
  @printf("B f2 ε = %.2e\n", abs(d4_f2(x) - B.d4_6(f2, x)))
  println()
  @printf("A f3 ε = %.2e\n", abs(d4_f3(x) - A.d4(f3, x)))
  @printf("B f3 ε = %.2e\n", abs(d4_f3(x) - B.d4_4(f3, x)))
  @printf("B f3 ε = %.2e\n", abs(d4_f3(x) - B.d4_6(f3, x)))
  println()
  @printf("A f4 ε = %.2e\n", abs(d4_f4(x) - A.d4(f4, x)))
  @printf("B f4 ε = %.2e\n", abs(d4_f4(x) - B.d4_4(f4, x)))
  @printf("B f4 ε = %.2e\n", abs(d4_f4(x) - B.d4_6(f4, x)))
  println()
  @printf("A f5 ε = %.2e\n", abs(d4_f5(x) - A.d4(f5, x)))
  @printf("B f5 ε = %.2e\n", abs(d4_f5(x) - B.d4_4(f5, x)))
  @printf("B f5 ε = %.2e\n", abs(d4_f5(x) - B.d4_6(f5, x)))
  println()
end

#=
# Todo: better performance tests
println("----- nested functions -----")
@time A.d4(f1, x)
@time A.d4(f2, x)
@time A.d4(f3, x)
@time A.d4(f4, x)
@time A.d4(f5, x)
println("----- direct formula (accuracy 4) -----")
@time B.d4_4(f1, x)
@time B.d4_4(f2, x)
@time B.d4_4(f3, x)
@time B.d4_4(f4, x)
@time B.d4_4(f5, x)
println("----- direct formula (accuracy 6) -----")
@time B.d4_6(f1, x)
@time B.d4_6(f2, x)
@time B.d4_6(f3, x)
@time B.d4_6(f4, x)
@time B.d4_6(f5, x)
=#

display([2, 4, 6, 8, 10] / 2)
