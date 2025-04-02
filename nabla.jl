include("deriv_rework.jl")

# ∇¹

function ∇1_2(f, x⃗, h::Float64=1e-3)
  h⃗ = [h for _ in x⃗]
  return (-(1 / 2)f(x⃗ - h⃗) + (1 / 2)f(x⃗ + h⃗)) / h
end

function ∇1_4(f, x⃗, h::Float64=1e-3)
  h⃗ = [h for _ in x⃗]
  return ((1 / 12)f(x⃗ - 2h⃗) - (2 / 3)f(x⃗ - h⃗) + (2 / 3)f(x⃗ + h⃗) - (1 / 12)f(x⃗ + 2h⃗)) / h
end

function ∇1_6(f, x⃗, h::Float64=1e-3)
  h⃗ = [h for _ in x⃗]
  return (-(1 / 60)f(x⃗ - 3h⃗) + (3 / 20)f(x⃗ - 2h⃗) - (3 / 4)f(x⃗ - h⃗) + (3 / 4)f(x⃗ + h⃗) - (3 / 20)f(x⃗ + 2h⃗) + (1 / 60)f(x⃗ + 3h⃗)) / h
end

function ∇1_8(f, x⃗, h::Float64=1e-3)
  h⃗ = [h for _ in x⃗]
  return ((1 / 280)f(x⃗ - 4h⃗) - (4 / 105)f(x⃗ - 3h⃗) + (1 / 5)f(x⃗ - 2h⃗) - (4 / 5)f(x⃗ - h⃗) + (4 / 5)f(x⃗ + h⃗) - (1 / 5)f(x⃗ + 2h⃗) + (4 / 105)f(x⃗ + 3h⃗) - (1 / 280)f(x⃗ + 4h⃗)) / h
end

# ∇²

function ∇2_2(f, x⃗, h::Float64=1e-3)
  h⃗ = [h for _ in x⃗]
  return (f(x⃗ - h⃗) - 2f(x⃗) + f(x⃗ + h⃗)) / h^2
end

function ∇2_4(f, x⃗, h::Float64=1e-3)
  h⃗ = [h for _ in x⃗]
  return (-(1 / 12)f(x⃗ - 2h⃗) + (4 / 3)f(x⃗ - h⃗) - (5 / 2)f(x⃗) + (4 / 3)f(x⃗ + h⃗) - (1 / 12)f(x⃗ + 2h⃗)) / h^2
end

function ∇2_6(f, x⃗, h::Float64=1e-3)
  h⃗ = [h for _ in x⃗]
  return ((1 / 90)f(x⃗ - 3h⃗) - (3 / 20)f(x⃗ - 2h⃗) + (3 / 2)f(x⃗ - h⃗) - (49 / 18)f(x⃗) + (3 / 2)f(x⃗ + h⃗) - (3 / 20)f(x⃗ + 2h⃗) + (1 / 90)f(x⃗ + 3h⃗)) / h^2
end

function ∇2_8(f, x⃗, h::Float64=1e-3)
  h⃗ = [h for _ in x⃗]
  return (-(1 / 560)f(x⃗ - 4h⃗) + (8 / 315)f(x⃗ - 3h⃗) - (1 / 5)f(x⃗ - 2h⃗) + (8 / 5)f(x⃗ - h⃗) - (205 / 72)f(x⃗) + (8 / 5)f(x⃗ + h⃗) - (1 / 5)f(x⃗ + 2h⃗) + (8 / 315)f(x⃗ + 3h⃗) - (1 / 560)f(x⃗ + 4h⃗)) / h^2
end

# ∇³

# ∇⁴

# ∇⁵

# ∇⁶
