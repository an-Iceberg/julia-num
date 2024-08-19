using Printf
include("integral.jl")
include("derivative.jl")

x::Float64 = -2
f(x) = ℯ^-(x - 2)^2 * sin(10x)
f⁽¹⁾(x) = -ℯ^-(x - 2)^2 * ((2x - 4)sin(10x) - 10cos(10x))

@printf "\033[4munderline\033[0m\n"

@printf "e^-(x-2)²⋅sin(10x)\n"
@printf "  d               = %.2e\n" abs(f⁽¹⁾(x) - d(1, f, x; max_prec=false))
@printf "  d max prec      = %.2e\n" abs(f⁽¹⁾(x) - d(1, f, x))
println()
for i in 0:10
  h = 1 / 10^i
  @printf "  d h=%.0e       = %.2e\n" h abs(f⁽¹⁾(x) - d(1, f, x; max_prec=false, h=h))
end
println()
@printf "  d h²            = %.2e\n" abs(f⁽¹⁾(x) - d(1, f, x; ext=true, max_prec=false))
@printf "  d h² h=0.1 d=5  = %.2e\n" abs(f⁽¹⁾(x) - d(1, f, x; h=0.1, depth=5, ext=true, max_prec=false))
@printf "  d h² h=0.01 d=3 = %.2e\n" abs(f⁽¹⁾(x) - d(1, f, x; h=0.01, depth=3, ext=true, max_prec=false))
@printf "  d h² max prec   = %.2e\n" abs(f⁽¹⁾(x) - d(1, f, x; ext=true))
println()

f₁(x) = x^2 - 2x + 1
f⁽¹⁾₁(x) = 2x - 2
@printf "x² - 2x + 1\n"
@printf "  d               = %.2e\n" abs(f⁽¹⁾₁(x) - d(1, f₁, x; max_prec=false))
@printf "  d max prec      = %.2e\n" abs(f⁽¹⁾₁(x) - d(1, f₁, x))
println()
for i in 0:10
  h = 1 / 10^i
  @printf "  d h=%.0e       = %.2e\n" h abs(f⁽¹⁾₁(x) - d(1, f₁, x; max_prec=false, h=h))
end
println()
@printf "  d h²            = %.2e\n" abs(f⁽¹⁾₁(x) - d(1, f₁, x; ext=true, max_prec=false))
@printf "  d h² max prec   = %.2e\n" abs(f⁽¹⁾₁(x) - d(1, f₁, x; ext=true))
println()

f₂(x) = 1 / x
f⁽¹⁾₂(x) = -1 / x^2
@printf "1/x\n"
@printf "  d               = %.2e\n" abs(f⁽¹⁾₂(x) - d(1, f₂, x; max_prec=false))
@printf "  d max prec      = %.2e\n" abs(f⁽¹⁾₂(x) - d(1, f₂, x))
println()
for i in 0:10
  h = 1 / 10^i
  @printf "  d h=%.0e       = %.2e\n" h abs(f⁽¹⁾₂(x) - d(1, f₂, x; max_prec=false, h=h))
end
println()
@printf "  d h²            = %.2e\n" abs(f⁽¹⁾₂(x) - d(1, f₂, x; ext=true, max_prec=false))
@printf "  d h² max prec   = %.2e\n" abs(f⁽¹⁾₂(x) - d(1, f₂, x; ext=true))
println()

f₃(x) = log(abs(x))
f⁽¹⁾₃(x) = 1 / x
@printf "log(|x|)\n"
@printf "  d               = %.2e\n" abs(f⁽¹⁾₃(x) - d(1, f₃, x; max_prec=false))
@printf "  d max prec      = %.2e\n" abs(f⁽¹⁾₃(x) - d(1, f₃, x))
println()
for i in 0:10
  h = 1 / 10^i
  @printf "  d h=%.0e       = %.2e\n" h abs(f⁽¹⁾₃(x) - d(1, f₃, x; max_prec=false, h=h))
end
println()
@printf "  d h²            = %.2e\n" abs(f⁽¹⁾₃(x) - d(1, f₃, x; ext=true, max_prec=false))
@printf "  d h² max prec   = %.2e\n" abs(f⁽¹⁾₃(x) - d(1, f₃, x; ext=true))
println()

f₄(x) = sin(x)
f⁽¹⁾₄(x) = cos(x)
@printf "sin(x)\n"
@printf "  d               = %.2e\n" abs(f⁽¹⁾₄(x) - d(1, f₄, x; max_prec=false))
@printf "  d max prec      = %.2e\n" abs(f⁽¹⁾₄(x) - d(1, f₄, x))
println()
for i in 0:10
  h = 1 / 10^i
  @printf "  d h=%.0e       = %.2e\n" h abs(f⁽¹⁾₄(x) - d(1, f₄, x; max_prec=false, h=h))
end
println()
@printf "  d h²            = %.2e\n" abs(f⁽¹⁾₄(x) - d(1, f₄, x; ext=true, max_prec=false))
@printf "  d h² max prec   = %.2e\n" abs(f⁽¹⁾₄(x) - d(1, f₄, x; ext=true))

# exit()
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

@printf "\033[4mpartial derivative\033[0m\n"
@printf "  ∂f⃗/∂x₁      ε = %.2e\n" abs(∂f⃗_∂x₁(x⃗[1]) - ∂(1, f⃗, x⃗, 1))
@printf "  ∂f⃗/∂x₁ h²   ε = %.2e\n" abs(∂f⃗_∂x₁(x⃗[1]) - ∂(1, f⃗, x⃗, 1; ext=true, max_prec=false))
@printf "  ∂f⃗/∂x₂      ε = %.2e\n" abs(∂f⃗_∂x₂(x⃗[2]) - ∂(1, f⃗, x⃗, 2))
@printf "  ∂f⃗/∂x₂ h²   ε = %.2e\n" abs(∂f⃗_∂x₂(x⃗[2]) - ∂(1, f⃗, x⃗, 2; ext=true, max_prec=false))
@printf "  ∂³f⃗/∂x₁³    ε = %.2e\n" abs(0 - ∂(3, f⃗, x⃗, 1))
@printf "  ∂³f⃗/∂x₁³ h² ε = %.2e\n" abs(0 - ∂(3, f⃗, x⃗, 1; ext=true, max_prec=false))
println()

@printf "\033[4mpartial derivative\033[0m\n"
@printf "  ∂f⃗/∂x₁      ε = %.2e\n" abs(∂f⃗_∂x₁(x⃗[1]) - ∂(1, f⃗, x⃗, 1))
@printf "  ∂f⃗/∂x₁ h²   ε = %.2e\n" abs(∂f⃗_∂x₁(x⃗[1]) - ∂(1, f⃗, x⃗, 1; ext=true, max_prec=false))
@printf "  ∂f⃗/∂x₂      ε = %.2e\n" abs(∂f⃗_∂x₂(x⃗[2]) - ∂(1, f⃗, x⃗, 2))
@printf "  ∂f⃗/∂x₂ h²   ε = %.2e\n" abs(∂f⃗_∂x₂(x⃗[2]) - ∂(1, f⃗, x⃗, 2; ext=true, max_prec=false))
@printf "  ∂²f⃗/∂x₂²    ε = %.2e\n" abs(∂²f⃗_∂x₂²(x⃗[2]) - ∂(2, f⃗, x⃗, 2))
@printf "  ∂²f⃗/∂x₂² h² ε = %.2e\n" abs(∂²f⃗_∂x₂²(x⃗[2]) - ∂(2, f⃗, x⃗, 2; ext=true, max_prec=false))
@printf "  ∂³f⃗/∂x₁³    ε = %.2e\n" abs(0 - ∂(3, f⃗, x⃗, 1))
@printf "  ∂³f⃗/∂x₁³ h² ε = %.2e\n" abs(0 - ∂(3, f⃗, x⃗, 1; ext=true, max_prec=false))
println()

println([round(x, digits=4) for x in ∇(1, f⃗, x⃗)])
println([round(x, digits=4) for x in ∇(2, f⃗, x⃗)])
println([round(x, digits=4) for x in ∇(3, f⃗, x⃗)])
println([round(x, digits=4) for x in ∇(4, f⃗, x⃗)])

f⃗₁(x⃗::Vector{<:Real})::Vector{<:Real} = [
  x⃗[1]^3 - x⃗[2]^2,
  x⃗[1]^2 - x⃗[2]
]
