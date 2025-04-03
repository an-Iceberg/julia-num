using Printf

include("deriv_rework.jl")
include("partial.jl")
include("nabla.jl")
include("int_rework.jl")

@printf "--------------------- Differentiation ---------------------\n\n"

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

@printf "--------------------- Partial differentiation ---------------------\n\n"

x⃗ = [Float64(π), Float64(π)]
f(x⃗) = sin(x⃗[1]) + cos(x⃗[2])
∂1x_f(x⃗) = cos(x⃗[1])
∂1y_f(x⃗) = -sin(x⃗[1])

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

@printf "--------------------- Gradient ---------------------\n\n"

x⃗ = [Float64(π), Float64(π)]
f(x⃗) = [sin(x⃗[1]), cos(x⃗[2])]
∂1x_f(x⃗) = cos(x⃗[1])
∂1y_f(x⃗) = -sin(x⃗[1])

@printf "∇1_2: %.2e\n" maximum([abs(analytical - numerical) for (analytical, numerical) in zip([∂1x_f(x⃗), ∂1y_f(x⃗)], ∇1_2(f, x⃗))])
