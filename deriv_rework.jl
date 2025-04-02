# https://en.wikipedia.org/wiki/Finite_difference_coefficient

# %%%%%%%%%%%%%%%%%%%% Derivative %%%%%%%%%%%%%%%%%%%%

# d¹
"""
Calculates the 1ˢᵗ order derivative of `f` at `x` with precision `h` using 2 points.
"""
function d1_2(f::Function, x::Real, h::Float64=1e-3)
  return (-0.5f(x - h) + 0.5f(x + h)) / h
end

function d1_4(f::Function, x::Real, h::Float64=1e-3)
  return ((1 / 12)f(x - 2h) - (2 / 3)f(x - h) + (2 / 3)f(x + h) - (1 / 12)f(x + 2h)) / h
end

function d1_6(f::Function, x::Real, h::Float64=1e-3)
  return (-(1 / 60)f(x - 3h) + 0.15f(x - 2h) - 0.75f(x - h) + 0.75f(x + h) - 0.15f(x + 2h) + (1 / 60)f(x + 3h)) / h
end

function d1_8(f::Function, x::Real, h::Float64=1e-3)
  return ((1 / 280)f(x - 4h) - (4 / 105)f(x - 3h) + 0.2f(x - 2h) - 0.8f(x - h) + 0.8f(x + h) - 0.2f(x + 2h) + (4 / 105)f(x + 3h) - (1 / 280)f(x + 4h)) / h
end

# d²

function d2_2(f::Function, x::Real, h::Float64=1e-3)
  return (f(x - h) - 2f(x) + f(x + h)) / h^2
end

function d2_4(f::Function, x::Real, h::Float64=1e-3)
  return (-(1 / 12)f(x - 2h) + (4 / 3)f(x - h) - (5 / 2)f(x) + (4 / 3)f(x + h) - (1 / 12)f(x + 2h)) / h^2
end

function d2_6(f::Function, x::Real, h::Float64=1e-3)
  return ((1 / 90)f(x - 3h) - (3 / 20)f(x - 2h) + (3 / 2)f(x - h) - (49 / 18)f(x) + (3 / 2)f(x + h) - (3 / 20)f(x + 2h) + (1 / 90)f(x + 3h)) / h^2
end

function d2_8(f::Function, x::Real, h::Float64=1e-3)
  return (-(1 / 560)f(x - 4h) + (8 / 315)f(x - 3h) - (1 / 5)f(x - 2h) + (8 / 5)f(x - h) - (205 / 72)f(x) + (8 / 5)f(x + h) - (1 / 5)f(x + 2h) + (8 / 315)f(x + 3h) - (1 / 560)f(x + 4h)) / h^2
end

# d³

function d3_2(f::Function, x::Real, h::Float64=1e-3)
  return (-(1 / 2)f(x - 2h) + f(x - h) - f(x + h) + (1 / 2)f(x + 2h)) / h^3
end

function d3_4(f::Function, x::Real, h::Float64=1e-3)
  return ((1 / 8)f(x - 3h) - f(x - 2h) + (13 / 8)f(x - h) - (13 / 8)f(x + h) + f(x + 2h) - (1 / 8)f(x + 3h)) / h^3
end

function d3_6(f::Function, x::Real, h::Float64=1e-3)
  return (-(7 / 240)f(x - 4h) + (3 / 10)f(x - 3h) - (169 / 120)f(x - 2h) + (61 / 30)f(x - h) - (61 / 30)f(x + h) + (169 / 120)f(x + 2h) - (3 / 10)f(x + 3h) + (7 / 240)f(x + 4h)) / h^3
end

# d⁴

function d4_2(f::Function, x::Real, h::Float64=1e-3)
  return (f(x - 2h) - 4f(x - h) + 6f(x) - 4f(x + h) + f(x + 2h)) / h^4
end

function d4_4(f::Function, x::Real, h::Float64=1e-3)
  return (-(1 / 6)f(x - 3h) + 2f(x - 2h) - (13 / 2)f(x - h) + (28 / 3)f(x) - (13 / 2)f(x + h) + 2f(x + 2h) - (1 / 6)f(x + 3h)) / h^4
end

function d4_6(f::Function, x::Real, h::Float64=1e-3)
  return ((7 / 240)f(x - 4h) - (2 / 5)f(x - 3h) + (169 / 60)f(x - 2h) - (122 / 15)f(x - h) + (91 / 8)f(x) - (122 / 15)f(x + h) + (169 / 60)f(x + 2h) - (2 / 5)f(x + 3h) + (7 / 240)f(x + 4h)) / h^4
end

# d⁵

function d5_2(f::Function, x::Real, h::Float64=1e-3)
  return (-(1 / 2)f(x - 3h) + 2f(x - 2h) - (5 / 2)f(x - h) + (5 / 2)f(x + h) - 2f(x + 2h) + (1 / 2)f(x + 3h)) / h^5
end

function d5_4(f::Function, x::Real, h::Float64=1e-3)
  return ((1 / 6)f(x - 4h) - (3 / 2)f(x - 3h) + (13 / 3)f(x - 2h) - (29 / 6)f(x - h) + (29 / 6)f(x + h) - (13 / 3)f(x + 2h) + (3 / 2)f(x + 3h) - (1 / 6)f(x + 4h)) / h^5
end

function d5_6(f::Function, x::Real, h::Float64=1e-3)
  return (-(13 / 288)f(x - 5h) + (19 / 36)f(x - 4h) - (87 / 32)f(x - 3h) + (13 / 2)f(x - 2h) - (323 / 48)f(x - h) + (323 / 48)f(x + h) - (13 / 2)f(x + 2h) + (87 / 32)f(x + 3h) - (19 / 36)f(x + 4h) + (13 / 288)f(x + 5h)) / h^5
end

# d⁶

function d6_2(f::Function, x::Real, h::Float64=1e-3)
  return (f(x - 3h) - 6f(x - 2h) + 15f(x - h) - 20f(x) + 15f(x + h) - 6f(x + 2h) + f(x + 3h)) / h^6
end

function d6_4(f::Function, x::Real, h::Float64=1e-3)
  return (-(1 / 4)f(x - 4h) + 3f(x - 3h) - 13f(x - 2h) + 29f(x - h) - (75 / 2)f(x) + 29f(x + h) - 13f(x + 2h) + 3f(x + 3h) - (1 / 4)f(x + 4h)) / h^6
end

function d6_6(f::Function, x::Real, h::Float64=1e-3)
  return ((13 / 240)f(x - 5h) - (19 / 24)f(x - 4h) + (87 / 16)f(x - 3h) - (39 / 2)f(x - 2h) + (323 / 8)f(x - h) - (1023 / 20)f(x) + (323 / 8)f(x - h) - (39 / 2)f(x - 2h) + (87 / 16)f(x - 3h) - (19 / 24)f(x - 4h) + (13 / 240)f(x - 5h)) / h^6
end
