# https://en.wikipedia.org/wiki/Finite_difference_coefficient

# https://discourse.julialang.org/t/how-to-do-partial-derivatives/19869/6
# https://stackoverflow.com/questions/54277219/partial-derivatives-in-julia
# https://discourse.julialang.org/t/round-a-float-and-force-a-specific-number-of-decimal-places/87547
# https://www.math.hkust.edu.hk/~mamu/courses/231/Slides/CH04_1B.pdf

# https://web.media.mit.edu/~crtaylor/calculator.html

# Rework idea: make use of this: d2(x) = d(x -> d(f, x), x)

# d¹
"""
Calculates the 1ˢᵗ order derivative df/dx of `f` at `x` with precision `h` using 2 points.
"""
function d1_2(f::Function, x::Real, h::Real=1e-3)::Real
  return (-f(x - h) + f(x + h)) / 2h
end

function d1_4(f::Function, x::Real, h::Real=1e-3)::Real
  return ((1 / 12)f(x - 2h) - (2 / 3)f(x - h) + (2 / 3)f(x + h) - (1 / 12)f(x + 2h)) / h
end

function d1_6(f::Function, x::Real, h::Real=1e-3)::Real
  return (
    -(1 / 60)f(x - 3h) + 0.15f(x - 2h) - 0.75f(x - h) + 0.75f(x + h) - 0.15f(x + 2h) +
    (1 / 60)f(x + 3h)
  ) / h
end

function d1_8(f::Function, x::Real, h::Real=1e-3)::Real
  return (
    (1 / 280)f(x - 4h) - (4 / 105)f(x - 3h) + 0.2f(x - 2h) - 0.8f(x - h) + 0.8f(x + h) -
    0.2f(x + 2h) + (4 / 105)f(x + 3h) - (1 / 280)f(x + 4h)
  ) / h
end

# d²

function d2_2(f::Function, x::Real, h::Real=1e-3)::Real
  return (f(x - h) - 2f(x) + f(x + h)) / h^2
end

function d2_4(f::Function, x::Real, h::Real=1e-3)::Real
  return (-(1 / 12)f(x - 2h) + (4 / 3)f(x - h) - 2.5f(x) + (4 / 3)f(x + h) - (1 / 12)f(x + 2h)) /
         h^2
end

function d2_6(f::Function, x::Real, h::Real=1e-3)::Real
  return (
    (1 / 90)f(x - 3h) - 0.15f(x - 2h) + 1.5f(x - h) - (49 / 18)f(x) + 1.5f(x + h) - 0.15f(x + 2h) +
    (1 / 90)f(x + 3h)
  ) / h^2
end

function d2_8(f::Function, x::Real, h::Real=1e-3)::Real
  return (
    -(1 / 560)f(x - 4h) + (8 / 315)f(x - 3h) - 0.2f(x - 2h) + 1.6f(x - h) - (205 / 72)f(x) +
    1.6f(x + h) - 0.2f(x + 2h) + (8 / 315)f(x + 3h) - (1 / 560)f(x + 4h)
  ) / h^2
end

# d³

function d3_2(f::Function, x::Real, h::Real=1e-3)::Real
  return (-0.5f(x - 2h) + f(x - h) - f(x + h) + 0.5f(x + 2h)) / h^3
end

function d3_4(f::Function, x::Real, h::Real=1e-3)::Real
  return (0.125f(x - 3h) - f(x - 2h) + 1.625f(x - h) - 1.625f(x + h) + f(x + 2h) - 0.125f(x + 3h)) /
         h^3
end

function d3_6(f::Function, x::Real, h::Real=1e-3)::Real
  return (
    -(7 / 240)f(x - 4h) + 0.3f(x - 3h) - (169 / 120)f(x - 2h) + (61 / 30)f(x - h) -
    (61 / 30)f(x + h) + (169 / 120)f(x + 2h) - 0.3f(x + 3h) + (7 / 240)f(x + 4h)
  ) / h^3
end

# d⁴

function d4_2(f::Function, x::Real, h::Real=1e-3)::Real
  return (f(x - 2h) - 4f(x - h) + 6f(x) - 4f(x + h) + f(x + 2h)) / h^4
end

function d4_4(f::Function, x::Real, h::Real=1e-3)::Real
  return (
    -(1 / 6)f(x - 3h) + 2f(x - 2h) - (13 / 2)f(x - h) + (28 / 3)f(x) - (13 / 2)f(x + h) +
    2f(x + 2h) - (1 / 6)f(x + 3h)
  ) / h^4
end

function d4_6(f::Function, x::Real, h::Real=1e-3)::Real
  return (
    (7 / 240)f(x - 4h) - 0.4f(x - 3h) + (169 / 60)f(x - 2h) - (122 / 15)f(x - h) + 11.375f(x) -
    (122 / 15)f(x + h) + (169 / 60)f(x + 2h) - 0.4f(x + 3h) + (7 / 240)f(x + 4h)
  ) / h^4
end

# d⁵

function d5_2(f::Function, x::Real, h::Real=1e-3)::Real
  return (-0.5f(x - 3h) + 2f(x - 2h) - 2.5f(x - h) + 2.5f(x + h) - 2f(x + 2h) + 0.5f(x + 3h)) / h^5
end

function d5_4(f::Function, x::Real, h::Real=1e-3)::Real
  return (
    (1 / 6)f(x - 4h) - (3 / 2)f(x - 3h) + (13 / 3)f(x - 2h) - (29 / 6)f(x - h) + (29 / 6)f(x + h) -
    (13 / 3)f(x + 2h) + (3 / 2)f(x + 3h) - (1 / 6)f(x + 4h)
  ) / h^5
end

function d5_6(f::Function, x::Real, h::Real=1e-3)::Real
  return (
    -(13 / 288)f(x - 5h) + (19 / 36)f(x - 4h) - 2.71875f(x - 3h) + 6.5f(x - 2h) -
    (323 / 48)f(x - h) + (323 / 48)f(x + h) - 6.5f(x + 2h) + 2.71875f(x + 3h) - (19 / 36)f(x + 4h) +
    (13 / 288)f(x + 5h)
  ) / h^5
end

# d⁶

function d6_2(f::Function, x::Real, h::Real=1e-3)::Real
  return (f(x - 3h) - 6f(x - 2h) + 15f(x - h) - 20f(x) + 15f(x + h) - 6f(x + 2h) + f(x + 3h)) / h^6
end

function d6_4(f::Function, x::Real, h::Real=1e-3)::Real
  return (
    -0.25f(x - 4h) + 3f(x - 3h) - 13f(x - 2h) + 29f(x - h) - 37.5f(x) + 29f(x + h) - 13f(x + 2h) +
    3f(x + 3h) - 0.25f(x + 4h)
  ) / h^6
end

# Todo: fix
function d6_6(f::Function, x::Real, h::Real=1e-3)::Real
  return (
    (13 / 240)f(x - 5h) - (19 / 24)f(x - 4h) + 5.4375f(x - 3h) - 19.5f(x - 2h) + 40.375f(x - h) -
    51.15f(x) + 40.375f(x + h) - 19.5f(x + 2h) + 5.4375f(x + 3h) - (19 / 24)f(x - 4h) +
    (13 / 240)f(x + 5h)
  ) / h^6
end
