using Printf
include("derivative.jl")
include("integral.jl")

x::Float64 = -2
f(x) = ℯ^-(x - 2)^2 * sin(10x)
f⁽¹⁾(x) = -ℯ^-(x - 2)^2 * ((2x - 4)sin(10x) - 10cos(10x))

@printf "\033[4munderline\033[0m\n"

@printf "e^-(x-2)²⋅sin(10x)\n"
@printf "  d h=1e-4      = %.2e\n" abs(f⁽¹⁾(x) - d(f, x; max_prec=false))
@printf "  d max prec    = %.2e\n" abs(f⁽¹⁾(x) - d(f, x))
@printf "  d h²          = %.2e\n" abs(f⁽¹⁾(x) - __d_h²(f, x))
@printf "  d h² max prec = %.2e\n" abs(f⁽¹⁾(x) - __d_h²_max_prec(f, x))
println()

f₁(x) = x^2 - 2x + 1
f⁽¹⁾₁(x) = 2x - 2
@printf "x² - 2x + 1\n"
@printf "  d h=1e-4      = %.2e\n" abs(f⁽¹⁾₁(x) - d(f₁, x; max_prec=false))
@printf "  d max prec    = %.2e\n" abs(f⁽¹⁾₁(x) - d(f₁, x))
@printf "  d h²          = %.2e\n" abs(f⁽¹⁾₁(x) - __d_h²(f₁, x))
@printf "  d h² max prec = %.2e\n" abs(f⁽¹⁾₁(x) - __d_h²_max_prec(f₁, x))
println()

f₂(x) = 1 / x
f⁽¹⁾₂(x) = -1 / x^2
@printf "1/x\n"
@printf "  d h=1e-4      = %.2e\n" abs(f⁽¹⁾₂(x) - d(f₂, x; max_prec=false))
@printf "  d max prec    = %.2e\n" abs(f⁽¹⁾₂(x) - d(f₂, x))
@printf "  d h²          = %.2e\n" abs(f⁽¹⁾₂(x) - __d_h²(f₂, x))
@printf "  d h² max prec = %.2e\n" abs(f⁽¹⁾₂(x) - __d_h²_max_prec(f₂, x))
println()

f₃(x) = log(abs(x))
f⁽¹⁾₃(x) = 1 / x
@printf "log(|x|)\n"
@printf "  d h=1e-4      = %.2e\n" abs(f⁽¹⁾₃(x) - d(f₃, x; max_prec=false))
@printf "  d max prec    = %.2e\n" abs(f⁽¹⁾₃(x) - d(f₃, x))
@printf "  d h²          = %.2e\n" abs(f⁽¹⁾₃(x) - __d_h²(f₃, x))
@printf "  d h² max prec = %.2e\n" abs(f⁽¹⁾₃(x) - __d_h²_max_prec(f₃, x))
println()

f₄(x) = sin(x)
f⁽¹⁾₄(x) = cos(x)
@printf "sin(x)\n"
@printf "  d h=1e-4      = %.2e\n" abs(f⁽¹⁾₄(x) - d(f₄, x; max_prec=false))
@printf "  d max prec    = %.2e\n" abs(f⁽¹⁾₄(x) - d(f₄, x))
@printf "  d h²          = %.2e\n" abs(f⁽¹⁾₄(x) - __d_h²(f₄, x))
@printf "  d h² max prec = %.2e\n" abs(f⁽¹⁾₄(x) - __d_h²_max_prec(f₄, x))
exit()
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
@printf "d midpoint    = %.2e\n" abs(d(f, x) - __d_midpoint(f, x))
@printf "d h² max prec = %.2e\n" abs(d(f, x) - __d_h²_max_prec(f, x))
# println("  --- depth = 2 ---")
# @printf "%.2e\n" abs(d(f, x) - __d_h²(f, x, 1e-1, 2))
# @printf "%.2e\n" abs(d(f, x) - __d_h²(f, x, 1e-2, 2))
# @printf "%.2e\n" abs(d(f, x) - __d_h²(f, x, 1e-3, 2))
# @printf "%.2e\n" abs(d(f, x) - __d_h²(f, x, 1e-4, 2))
# println("  --- depth = 3 ---")
# @printf "%.2e\n" abs(d(f, x) - __d_h²(f, x, 1e-1, 3))
# @printf "%.2e\n" abs(d(f, x) - __d_h²(f, x, 1e-2, 3))
# @printf "%.2e\n" abs(d(f, x) - __d_h²(f, x, 1e-3, 3))
# @printf "%.2e\n" abs(d(f, x) - __d_h²(f, x, 1e-4, 3))
# println("  --- depth = 4 ---")
# @printf "%.2e\n" abs(d(f, x) - __d_h²(f, x, 1e-1, 4))
# @printf "%.2e\n" abs(d(f, x) - __d_h²(f, x, 1e-2, 4))
# @printf "%.2e\n" abs(d(f, x) - __d_h²(f, x, 1e-3, 4))
# @printf "%.2e\n" abs(d(f, x) - __d_h²(f, x, 1e-4, 4))
# println("  --- depth = 5 ---")
# @printf "%.2e\n" abs(d(f, x) - __d_h²(f, x, 1e-1, 5))
# @printf "%.2e\n" abs(d(f, x) - __d_h²(f, x, 1e-2, 5))
# @printf "%.2e\n" abs(d(f, x) - __d_h²(f, x, 1e-3, 5))
# @printf "%.2e\n" abs(d(f, x) - __d_h²(f, x, 1e-4, 5))

@printf "d³ midpoint    = %.2e\n" abs(d³(f, x) - __d³_midpoint(f, x))
@printf "d³ h²          = %.2e\n" abs(d³(f, x) - __d³_h²(f, x))
@printf "d³ h² max prec = %.2e\n" abs(d³(f, x) - __d³_h²_max_prec(f, x))
@printf "d⁴ midpoint    = %.2e\n" abs(d⁴(f, x) - __d⁴_midpoint(f, x))
@printf "d⁴ h²          = %.2e\n" abs(d⁴(f, x) - __d⁴_h²(f, x))
@printf "d⁴ h² max prec = %.2e\n" abs(d⁴(f, x) - __d⁴_h²_max_prec(f, x))

