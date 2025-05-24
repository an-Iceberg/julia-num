include("src/gradient.jl")
include("src/jacobian.jl")

f1(x) = x[1] + x[2]^2 - x[3]^2 - 13
f2(x) = log(x[2] / 4) + exp(0.5x[3] - 1) - 1
f3(x) = (x[2] - 3)^2 - x[3]^3 + 7

f(λ) = [
  f1(λ),
  f2(λ),
  f3(λ),
]
x⃗ = [1.5, 3, 2.5]

print("∇: ")
display(∇(f, x⃗))
print("D: ")
display(D1_4(f, x⃗))
