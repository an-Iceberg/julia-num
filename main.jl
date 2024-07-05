using Printf
include("derivative.jl")
include("integral.jl")

x::Float64 = 2
f(x) = ℯ^-(x - 2)^2 * sin(10x)
f⁽¹⁾(x) = -ℯ^-(x - 2)^2 * ((2x - 4)sin(10x) - 10cos(10x))

@printf "\033[4mapproximation errors\033[0m\n"
@printf "  5 point          = %.2e\n" abs(f⁽¹⁾(x) - d(f, x, max_prec=false))
@printf "  h²               = %.2e\n" abs(f⁽¹⁾(x) - __d_h²(f, x))
@printf "  5 point max prec = %.2e\n" abs(f⁽¹⁾(x) - d(f, x))

# Todo: performance testing

# https://discourse.julialang.org/t/how-to-do-partial-derivatives/19869/6
# https://stackoverflow.com/questions/54277219/partial-derivatives-in-julia
# https://discourse.julialang.org/t/round-a-float-and-force-a-specific-number-of-decimal-places/87547
# https://www.math.hkust.edu.hk/~mamu/courses/231/Slides/CH04_1B.pdf

println()

x⃗::Vector{Real} = [2, 2]

f⃗(x⃗::Vector{<:Real})::Real = 0.5x⃗[1]^2 + x⃗[2]^3
∂f⃗_∂x₁(x::Real)::Real = x
∂f⃗_∂x₂(x::Real)::Real = 3x^2
∂²f⃗_∂x₂²(x::Real)::Real = 6x

@printf "\033[4mapproximation errors\033[0m\n"
@printf "  5 point ∂f⃗/∂x₁   = %.2e\n" abs(∂f⃗_∂x₁(x⃗[1]) - ∂(f⃗, x⃗, 1, max_prec=false))
@printf "  5 point ∂f⃗/∂x₂   = %.2e\n" abs(∂f⃗_∂x₂(x⃗[2]) - ∂(f⃗, x⃗, 2, max_prec=false))
@printf "  5 point ∂³f⃗/∂x₁³ = %.2e\n" abs(0 - ∂³(f⃗, x⃗, 1, max_prec=false))
println()
@printf "\033[4mmax precision\033[0m\n"
@printf "  5 point ∂f⃗/∂x₁   = %.2e\n" abs(∂f⃗_∂x₁(x⃗[1]) - ∂(f⃗, x⃗, 1))
@printf "  5 point ∂f⃗/∂x₂   = %.2e\n" abs(∂f⃗_∂x₂(x⃗[2]) - ∂(f⃗, x⃗, 2))
@printf "  5 point ∂²f⃗/∂x₂² = %.2e\n" abs(∂²f⃗_∂x₂²(x⃗[2]) - ∂²(f⃗, x⃗, 2))
@printf "  5 point ∂³f⃗/∂x₁³ = %.2e\n" abs(0 - ∂³(f⃗, x⃗, 1))

println()

println([round(x, digits=4) for x in ∇(f⃗, x⃗)])
println([round(x, digits=4) for x in ∇²(f⃗, x⃗)])
println([round(x, digits=4) for x in ∇³(f⃗, x⃗)])
println([round(x, digits=4) for x in ∇⁴(f⃗, x⃗)])

println()
@printf "%.2e\n" abs(d(f, x) - __d_midpoint(f, x))
@printf "%.2e\n" abs(d³(f, x) - __d³_midpoint(f, x))
@printf "%.2e\n" abs(d³(f, x) - __d³_h²(f, x))
@printf "%.2e\n" abs(d⁴(f, x) - __d⁴_midpoint(f, x))
@printf "%.2e\n" abs(d⁴(f, x) - __d⁴_h²(f, x))
