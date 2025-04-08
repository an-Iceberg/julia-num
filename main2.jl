using Printf

include("derivative.jl")
include("partial.jl")
include("gradient.jl")
include("jacobian.jl")
include("integral.jl")

@printf "--------------------- Differentiation (df) ---------------------\n\n"

x = Float64(π)
f(x) = sin(x)
d1f(x) = cos(x)
d2f(x) = -sin(x)
d3f(x) = -cos(x)
d4f(x) = sin(x)
d5f(x) = cos(x)
d6f(x) = -sin(x)

@printf "d1_2: %.2e\n" abs(d1f(x) - d1_2(f, x))
@printf "d1_4: %.2e\n" abs(d1f(x) - d1_4(f, x))
@printf "d1_6: %.2e\n" abs(d1f(x) - d1_6(f, x))
@printf "d1_8: %.2e\n" abs(d1f(x) - d1_8(f, x))
@printf "\n"
@printf "d2_2: %.2e\n" abs(d2f(x) - d2_2(f, x))
@printf "d2_4: %.2e\n" abs(d2f(x) - d2_4(f, x))
@printf "d2_6: %.2e\n" abs(d2f(x) - d2_6(f, x))
@printf "d2_8: %.2e\n" abs(d2f(x) - d2_8(f, x))
@printf "\n"
@printf "d3_2: %.2e\n" abs(d3f(x) - d3_2(f, x))
@printf "d3_4: %.2e\n" abs(d3f(x) - d3_4(f, x))
@printf "d3_6: %.2e\n" abs(d3f(x) - d3_6(f, x))
@printf "\n"
@printf "d4_2: %.2e\n" abs(d4f(x) - d4_2(f, x, 0.01))
@printf "d4_4: %.2e\n" abs(d4f(x) - d4_4(f, x, 0.01))
@printf "d4_6: %.2e\n" abs(d4f(x) - d4_6(f, x, 0.01))
@printf "\n"
@printf "d5_2: %.2e\n" abs(d5f(x) - d5_2(f, x, 0.1))
@printf "d5_4: %.2e\n" abs(d5f(x) - d5_4(f, x, 0.1))
@printf "d5_6: %.2e\n" abs(d5f(x) - d5_6(f, x, 0.1))
@printf "\n"
@printf "d6_2: %.2e\n" abs(d6f(x) - d6_2(f, x, 1.0))
@printf "d6_4: %.2e\n" abs(d6f(x) - d6_4(f, x, 1.0))
@printf "d6_6: %.2e\n" abs(d6f(x) - d6_6(f, x, 1.0))
@printf "\n"

@printf "--------------------- Partial differentiation (∂f) ---------------------\n\n"

x⃗ = [Float64(π), Float64(π)]
f(x⃗) = sin(x⃗[1]) + cos(x⃗[2])
∂1x_f(x⃗) = cos(x⃗[1])
∂1y_f(x⃗) = -sin(x⃗[1])
∂2x_f(x⃗) = -sin(x⃗[1])
∂2y_f(x⃗) = -cos(x⃗[1])
∂3x_f(x⃗) = -cos(x⃗[1])
∂3y_f(x⃗) = sin(x⃗[1])
∂4x_f(x⃗) = sin(x⃗[1])
∂4y_f(x⃗) = cos(x⃗[1])
∂5x_f(x⃗) = cos(x⃗[1])
∂5y_f(x⃗) = -sin(x⃗[1])
∂6x_f(x⃗) = -sin(x⃗[1])
∂6y_f(x⃗) = -cos(x⃗[1])

println(∂1_2(f, x⃗, 1))
println()

@printf "∂1_2: %.2e\n" abs(∂1x_f(x⃗) - ∂1_2(f, x⃗, 1))
@printf "∂1_4: %.2e\n" abs(∂1x_f(x⃗) - ∂1_4(f, x⃗, 1))
@printf "∂1_6: %.2e\n" abs(∂1x_f(x⃗) - ∂1_6(f, x⃗, 1))
@printf "∂1_8: %.2e\n" abs(∂1x_f(x⃗) - ∂1_8(f, x⃗, 1))
@printf "\n"
@printf "∂1_2: %.2e\n" abs(∂1y_f(x⃗) - ∂1_2(f, x⃗, 2))
@printf "∂1_4: %.2e\n" abs(∂1y_f(x⃗) - ∂1_4(f, x⃗, 2))
@printf "∂1_6: %.2e\n" abs(∂1y_f(x⃗) - ∂1_6(f, x⃗, 2))
@printf "∂1_8: %.2e\n" abs(∂1y_f(x⃗) - ∂1_8(f, x⃗, 2))
@printf "\n"
@printf "∂2_2: %.2e\n" abs(∂2x_f(x⃗) - ∂2_2(f, x⃗, 1))
@printf "∂2_4: %.2e\n" abs(∂2x_f(x⃗) - ∂2_4(f, x⃗, 1))
@printf "∂2_6: %.2e\n" abs(∂2x_f(x⃗) - ∂2_6(f, x⃗, 1))
@printf "∂2_8: %.2e\n" abs(∂2x_f(x⃗) - ∂2_8(f, x⃗, 1))
@printf "\n"
@printf "∂2_2: %.2e\n" abs(∂2y_f(x⃗) - ∂2_2(f, x⃗, 2))
@printf "∂2_4: %.2e\n" abs(∂2y_f(x⃗) - ∂2_4(f, x⃗, 2))
@printf "∂2_6: %.2e\n" abs(∂2y_f(x⃗) - ∂2_6(f, x⃗, 2))
@printf "∂2_8: %.2e\n" abs(∂2y_f(x⃗) - ∂2_8(f, x⃗, 2))
@printf "\n"
@printf "∂3_2: %.2e\n" abs(∂3x_f(x⃗) - ∂3_2(f, x⃗, 1))
@printf "∂3_4: %.2e\n" abs(∂3x_f(x⃗) - ∂3_4(f, x⃗, 1))
@printf "∂3_6: %.2e\n" abs(∂3x_f(x⃗) - ∂3_6(f, x⃗, 1))
@printf "\n"
@printf "∂3_2: %.2e\n" abs(∂3y_f(x⃗) - ∂3_2(f, x⃗, 2))
@printf "∂3_4: %.2e\n" abs(∂3y_f(x⃗) - ∂3_4(f, x⃗, 2))
@printf "∂3_6: %.2e\n" abs(∂3y_f(x⃗) - ∂3_6(f, x⃗, 2))
@printf "\n"
@printf "∂4_2: %.2e\n" abs(∂4x_f(x⃗) - ∂4_2(f, x⃗, 1, 0.1))
@printf "∂4_4: %.2e\n" abs(∂4x_f(x⃗) - ∂4_4(f, x⃗, 1, 0.1))
@printf "∂4_6: %.2e\n" abs(∂4x_f(x⃗) - ∂4_6(f, x⃗, 1, 0.1))
@printf "\n"
@printf "∂4_2: %.2e\n" abs(∂4y_f(x⃗) - ∂4_2(f, x⃗, 2, 0.1))
@printf "∂4_4: %.2e\n" abs(∂4y_f(x⃗) - ∂4_4(f, x⃗, 2, 0.1))
@printf "∂4_6: %.2e\n" abs(∂4y_f(x⃗) - ∂4_6(f, x⃗, 2, 0.1))
@printf "\n"
@printf "∂5_2: %.2e\n" abs(∂5x_f(x⃗) - ∂5_2(f, x⃗, 1, 0.1))
@printf "∂5_4: %.2e\n" abs(∂5x_f(x⃗) - ∂5_4(f, x⃗, 1, 0.1))
@printf "∂5_6: %.2e\n" abs(∂5x_f(x⃗) - ∂5_6(f, x⃗, 1, 0.1))
@printf "\n"
@printf "∂5_2: %.2e\n" abs(∂5y_f(x⃗) - ∂5_2(f, x⃗, 2, 0.1))
@printf "∂5_4: %.2e\n" abs(∂5y_f(x⃗) - ∂5_4(f, x⃗, 2, 0.1))
@printf "∂5_6: %.2e\n" abs(∂5y_f(x⃗) - ∂5_6(f, x⃗, 2, 0.1))
@printf "\n"
@printf "∂6_2: %.2e\n" abs(∂6x_f(x⃗) - ∂6_2(f, x⃗, 1, 0.1))
@printf "∂6_4: %.2e\n" abs(∂6x_f(x⃗) - ∂6_4(f, x⃗, 1, 0.1))
@printf "∂6_6: %.2e\n" abs(∂6x_f(x⃗) - ∂6_6(f, x⃗, 1, 0.1))
@printf "\n"
@printf "∂6_2: %.2e\n" abs(∂6y_f(x⃗) - ∂6_2(f, x⃗, 2, 0.1))
@printf "∂6_4: %.2e\n" abs(∂6y_f(x⃗) - ∂6_4(f, x⃗, 2, 0.1))
@printf "∂6_6: %.2e\n" abs(∂6y_f(x⃗) - ∂6_6(f, x⃗, 2, 0.1))
@printf "\n"

@printf "--------------------- Gradient (∇f) ---------------------\n\n"

x⃗ = [Float64(π), Float64(π)]
f(x⃗) = [sin(x⃗[1]), cos(x⃗[2])]
∂1x_f(x⃗) = cos(x⃗[1])
∂1y_f(x⃗) = -sin(x⃗[1])
∂2x_f(x⃗) = -sin(x⃗[1])
∂2y_f(x⃗) = -cos(x⃗[1])
∂3x_f(x⃗) = -cos(x⃗[1])
∂3y_f(x⃗) = sin(x⃗[1])
∂4x_f(x⃗) = sin(x⃗[1])
∂4y_f(x⃗) = cos(x⃗[1])
∂5x_f(x⃗) = cos(x⃗[1])
∂5y_f(x⃗) = -sin(x⃗[1])
∂6x_f(x⃗) = -sin(x⃗[1])
∂6y_f(x⃗) = -cos(x⃗[1])

# a = analytical solution, n = numerical solution
@printf "∇1_2: %.2e\n" maximum([abs(a - n) for (a, n) in zip([∂1x_f(x⃗), ∂1y_f(x⃗)], ∇1_2(f, x⃗))])
@printf "∇1_4: %.2e\n" maximum([abs(a - n) for (a, n) in zip([∂1x_f(x⃗), ∂1y_f(x⃗)], ∇1_4(f, x⃗))])
@printf "∇1_6: %.2e\n" maximum([abs(a - n) for (a, n) in zip([∂1x_f(x⃗), ∂1y_f(x⃗)], ∇1_6(f, x⃗))])
@printf "∇1_8: %.2e\n" maximum([abs(a - n) for (a, n) in zip([∂1x_f(x⃗), ∂1y_f(x⃗)], ∇1_8(f, x⃗))])
@printf "\n"
@printf "∇2_2: %.2e\n" maximum([abs(a - n) for (a, n) in zip([∂2x_f(x⃗), ∂2y_f(x⃗)], ∇2_2(f, x⃗))])
@printf "∇2_4: %.2e\n" maximum([abs(a - n) for (a, n) in zip([∂2x_f(x⃗), ∂2y_f(x⃗)], ∇2_4(f, x⃗))])
@printf "∇2_6: %.2e\n" maximum([abs(a - n) for (a, n) in zip([∂2x_f(x⃗), ∂2y_f(x⃗)], ∇2_6(f, x⃗))])
@printf "∇2_8: %.2e\n" maximum([abs(a - n) for (a, n) in zip([∂2x_f(x⃗), ∂2y_f(x⃗)], ∇2_8(f, x⃗))])
@printf "\n"
@printf "∇3_2: %.2e\n" maximum([abs(a - n) for (a, n) in zip([∂3x_f(x⃗), ∂3y_f(x⃗)], ∇3_2(f, x⃗))])
@printf "∇3_4: %.2e\n" maximum([abs(a - n) for (a, n) in zip([∂3x_f(x⃗), ∂3y_f(x⃗)], ∇3_4(f, x⃗))])
@printf "∇3_6: %.2e\n" maximum([abs(a - n) for (a, n) in zip([∂3x_f(x⃗), ∂3y_f(x⃗)], ∇3_6(f, x⃗))])
@printf "\n"
@printf "∇4_2: %.2e\n" maximum([abs(a - n) for (a, n) in zip([∂4x_f(x⃗), ∂4y_f(x⃗)], ∇4_2(f, x⃗))])
@printf "∇4_4: %.2e\n" maximum([abs(a - n) for (a, n) in zip([∂4x_f(x⃗), ∂4y_f(x⃗)], ∇4_4(f, x⃗))])
@printf "∇4_6: %.2e\n" maximum([abs(a - n) for (a, n) in zip([∂4x_f(x⃗), ∂4y_f(x⃗)], ∇4_6(f, x⃗))])
@printf "\n"
@printf "∇5_2: %.2e\n" maximum([abs(a - n) for (a, n) in zip([∂5x_f(x⃗), ∂5y_f(x⃗)], ∇5_2(f, x⃗, 0.01))])
@printf "∇5_4: %.2e\n" maximum([abs(a - n) for (a, n) in zip([∂5x_f(x⃗), ∂5y_f(x⃗)], ∇5_4(f, x⃗, 0.01))])
@printf "∇5_6: %.2e\n" maximum([abs(a - n) for (a, n) in zip([∂5x_f(x⃗), ∂5y_f(x⃗)], ∇5_6(f, x⃗, 0.01))])
@printf "\n"
@printf "∇6_2: %.2e\n" maximum([abs(a - n) for (a, n) in zip([∂6x_f(x⃗), ∂6y_f(x⃗)], ∇6_2(f, x⃗, 0.01))])
@printf "∇6_4: %.2e\n" maximum([abs(a - n) for (a, n) in zip([∂6x_f(x⃗), ∂6y_f(x⃗)], ∇6_4(f, x⃗, 0.01))])
@printf "∇6_6: %.2e\n" maximum([abs(a - n) for (a, n) in zip([∂6x_f(x⃗), ∂6y_f(x⃗)], ∇6_6(f, x⃗, 0.01))])

@printf "\n"

@printf "--------------------- Jacobian (Df) ---------------------\n\n"

f1(x) = x[1] + x[2]^2 - x[3]^2 - 13
f2(x) = log(x[2] / 4) + ℯ^(0.5x[3] - 1) - 1
f3(x) = (x[2] - 3)^2 - x[3]^3 + 7

f(λ) = [f1(λ), f2(λ), f3(λ)]
x⃗ = [1.5, 3, 2.5]

display(D1_2(f, x⃗))
display(D1_4(f, x⃗))
display(D1_6(f, x⃗))
display(D1_8(f, x⃗))

@printf "--------------------- Integral (∫f) ---------------------\n\n"

f(x) = cos(x)
F(x) = sin(x)
a = 0.0
b = 5.0

@printf "Simpson's ⅓:  ε = %.2e\n" abs(∫_2(a, b, f) - (F(b) - F(a)))
@printf "Simpson's ⅜:  ε = %.2e\n" abs(∫_3(a, b, f) - (F(b) - F(a)))
@printf "Boole's :     ε = %.2e\n" abs(∫_4(a, b, f) - (F(b) - F(a)))
