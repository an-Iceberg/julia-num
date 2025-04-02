using Printf

# %%%%%%%%%%%%%%%%%%%% Integration %%%%%%%%%%%%%%%%%%%%

@enum IntMethod simpsons_1_3 simpsons_3_8 booles

function simpson_⅓(a::Real, b::Real, f::Function, h::Real=0.001)::Real
  n = ceil((b - a) / h)
  x(i) = a + (i * h)

  Σ₁ = sum([f(x(2i - 1)) for i in 1:(n / 2)])
  Σ₂ = sum([f(x(2i)) for i in 1:((n / 2) - 1)])

  return (1 / 3) * h * (f(a) + 4Σ₁ + 2Σ₂ + f(b))
end

function simpson_⅜(a::Real, b::Real, f::Function, h::Real=0.001)::Real
  n = ceil((b - a) / h)
  x(i) = a + (i * h)

  Σ₁ = collect(1:n)
  filter!(i -> i % 3 != 0, Σ₁)
  Σ₁ = sum([f(x(i)) for i in Σ₁])
  Σ₂ = sum([f(x(3i)) for i in 1:((n / 3) - 1)])

  return (3 / 8) * h * (f(a) + 3Σ₁ + 2Σ₂ + f(b))

  # Σ = f(a) + f(b)
  # for i in 1:n
  #   if i % 3 == 0
  #     Σ += 2 * f(x(i))
  #   else
  #     Σ += 3 * f(x(i))
  #   end
  # end
  # return (3 / 8) * h * Σ

  # return (3 / 8) * h *
  #        sum([f(x(3i - 3)) + f(x(3i - 2)) + f(x(3i - 1)) + f(x(3i)) for i in 1:(n / 3 - 1)])
end

function boole(a::Real, b::Real, f::Function, h::Real=0.001)::Real
  n = ceil((b - a) / h)
  x(i) = a + (i * h)

  Σ₁ = sum([f(x(1 + 2i)) for i in 0:((n / 2) - 1)])
  Σ₂ = sum([f(x(2 + 4i)) for i in 0:((n / 4) - 1)])
  Σ₃ = sum([f(x(4 + 4i)) for i in 0:((n / 4) - 2)])

  return (2 / 45) * h * (7 * (f(a) + f(b)) + 32Σ₁ + 12Σ₂ + 14Σ₃)
end

function ∫(a::Real, b::Real, f::Function, h::Real=0.001, method::IntMethod=boole)::Real
  if method == simpsons_1_3
    return simpson_⅓(a, b, f, h)
  elseif method == simpsons_3_8
    return simpson_⅜(a, b, f, h)
  else
    return boole(a, b, f, h)
  end
end

f(x) = cos(x)
F(x) = sin(x)
a = 0.0
b = 5.0
h = 1e-2

@printf "Simpson's ⅓:  ε = %.2e  h = %.0e\n" abs(∫(a, b, f, h, simpsons_1_3) - (F(b) - F(a))) h
@printf "Simpson's ⅜:  ε = %.2e  h = %.0e\n" abs(∫(a, b, f, h, simpsons_3_8) - (F(b) - F(a))) h
@printf "Boole's :     ε = %.2e  h = %.0e\n" abs(∫(a, b, f, h, booles) - (F(b) - F(a))) h

# %%%%%%%%%%%%%%%%%%%% Differentiation %%%%%%%%%%%%%%%%%%%%

function d(degree::UInt, f::Function, x::Real, h::Real=0.001, acc::UInt=6)::Real
end

# %%%%%%%%%%%%%%%%%%%% Partial Differentiation %%%%%%%%%%%%%%%%%%%%

function ∂(degree::UInt, f::Function, x::Vector{<:Real}, i::UInt, h::Real=0.001, acc::UInt=6)::Real
end

# %%%%%%%%%%%%%%%%%%%% Partial Differentiation %%%%%%%%%%%%%%%%%%%%

function ∇()
end
