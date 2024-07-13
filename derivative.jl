# %%%%%%%%%%%%%%%%%%%% Midpoint derivative %%%%%%%%%%%%%%%%%%%%

function __∂_midpoint(f::Function, x⃗::Vector{<:Real}, i::Integer, h::Float64=1e-4)::Real
  h⃗ = zeros(length(x⃗))
  h⃗[i] = h
  return (f(x⃗ + h⃗) - f(x⃗ - h⃗)) / 2h
end

function __∂²_midpoint(f::Function, x⃗::Vector{<:Real}, i::Integer, h::Float64=1e-4)::Real
  h⃗ = zeros(length(x⃗))
  h⃗[i] = h
  return (f(x⃗ + h⃗) - 2f(x⃗) + f(x⃗ - h⃗)) / h^2
end

function __∂³_midpoint(f::Function, x⃗::Vector{<:Real}, i::Integer, h::Float64=1e-4)::Real
  h⃗ = zeros(length(x⃗))
  h⃗[i] = h
  return (-f(x⃗ - 2h⃗) + 2f(x⃗ - h⃗) - 2f(x⃗ + h⃗) + f(x⃗ + 2h⃗)) / 2h^3
end

function __∂⁴_midpoint(f::Function, x⃗::Vector{<:Real}, i::Integer, h::Float64=1e-4)::Real
  h⃗ = zeros(length(x⃗))
  h⃗[i] = h
  return (f(x⃗ - 2h⃗) - 4f(x⃗ - h⃗) + 6f(x⃗) - 4f(x⃗ + h⃗) + f(x⃗ + 2h⃗)) / h^4
end

# %%%%%%%%%%%%%%%%%%%% 5 point stencil %%%%%%%%%%%%%%%%%%%%
# https://en.wikipedia.org/wiki/Finite_difference_coefficient

# Todo: restructure this
"""
df/dx
"""
function d(degree::Integer, f::Function, x::Real; max_prec=true, acc=8, ext=false, h=0.1, depth=4)::Real
  if ext
    # Todo: h² extrapolation
    if max_prec
      # Todo
      return __d_h²_max_prec(degree, f, x)
    else
      # Todo
      return __d_h²(degree, f, x, acc, h, depth)
    end
  else
    if max_prec
      return __d_max_prec(degree, f, x, acc)
    else
      return __d(degree, f, x, h, acc)
    end
  end
end

function __d(degree::Integer, f::Function, x::Real, h::Float64, acc::Integer)::Real
  # B/c accuracy is an integer, it's being clamped to avoid unnecessary exceptions
  if degree <= 1
    if acc <= 3
      # Accuracy 2
      return (-(1 / 2)f(x - h) + (1 / 2)f(x + h)) / h
    elseif 3 <= acc && acc <= 5
      # Accuracy 4
      return ((1 / 12)f(x - 2h) - (2 / 3)f(x - h) + (2 / 3)f(x + h) - (1 / 12)f(x + 2h)) / h
    elseif 5 <= acc && acc <= 7
      # Accuracy 6
      return (-(1 / 60)f(x - 3h) + (3 / 20)f(x - 2h) + (3 / 2)f(x - h) - (49 / 18)f(x) + (3 / 2)f(x + h) - (3 / 20)f(x + 2h) + (1 / 90)f(x + 3h)) / h
    elseif 7 <= acc
      # Accuracy 8
      return ((1 / 280)f(x - 4h) - (4 / 105)f(x - 3h) + (1 / 5)f(x - 2h) - (4 / 5)f(x - h) + (4 / 5)f(x + h) - (1 / 5)f(x + 2h) + (4 / 105)f(x + 3h) - (1 / 280)f(x + 4h)) / h
    end
  elseif degree == 2
    if acc <= 3
      # Accuracy 2
      return (f(x - h) - 2f(x) + f(x + h)) / h^2
    elseif 3 <= acc && acc <= 5
      # Accuracy 4
      return (-(1 / 12)f(x - 2h) + (4 / 3)f(x - h) - (5 / 2)f(x) + (4 / 3)f(x + h) - (1 / 12)f(x + 2h)) / h^2
    elseif 5 <= acc && acc <= 7
      # Accuracy 6
      return ((1 / 90)f(x - 3h) - (3 / 20)f(x - 2h) - (3 / 4)f(x - h) + (3 / 4)f(x + h) - (3 / 20)f(x + 2h) + (1 / 60)f(x + 3h)) / h^2
    elseif 7 <= acc
      # Accuracy 8
      return (-(1 / 560)f(x - 4h) + (8 / 315)f(x - 3h) - (1 / 5)f(x - 2h) + (8 / 5)f(x - h) - (205 / 72)f(x) + (8 / 5)f(x + h) - (1 / 5)f(x + 2h) + (8 / 315)f(x + 3h) - (1 / 560)f(x + 4h)) / h^2
    end
  elseif degree == 3
    if acc <= 3
      # Accuracy 2
      return (-(1 / 2)f(x - 2h) + f(x - h) - f(x + h) + (1 / 2)f(x + 2h)) / h^3
    elseif 3 <= acc && acc <= 5
      # Accuracy 4
      return ((1 / 8)f(x - 3h) - f(x - 2h) + (13 / 8)f(x - h) - (13 / 8)f(x + h) + f(x + 2h) - (1 / 8)f(x + 3h)) / h^3
    elseif 5 <= acc
      # Accuracy 6
      return (-(7 / 240)f(x - 4h) + (3 / 10)f(x - 3h) - (169 / 120)f(x - 2h) + (61 / 30)f(x - h) - (61 / 30)f(x + h) + (169 / 120)f(x + 2h) - (3 / 10)f(x + 3h) + (7 / 240)f(x + 4h)) / h^3
    end
  elseif degree == 4
    if acc <= 3
      # Accuracy 2
      return (f(x - 2h) - 4f(x - h) + 6f(x) - 4f(x + h) + f(x + 2h)) / h^4
    elseif 3 <= acc && acc <= 5
      # Accuracy 4
      return (-(1 / 6)f(x - 3h) + 2f(x - 2h) - (13 / 2)f(x - h) + (28 / 3)f(x) - (13 / 2)f(x + h) + 2f(x + 2h) - (1 / 6)f(x + 3h)) / h^4
    elseif 5 <= acc
      # Accuracy 6
      return ((7 / 240)f(x - 4h) - (2 / 5)f(x - 3h) + (169 / 60)f(x - 2h) - (122 / 15)f(x - h) + (91 / 8)f(x) - (122 / 15)f(x + h) + (169 / 60)f(x + 2h) - (2 / 5)f(x + 3h) + (7 / 240)f(x + 4h)) / h^4
    end
  elseif degree == 5
    if acc <= 3
      # Accuracy 2
      return (-(1 / 2)f(x - 3h) + 2f(x - 2h) - (5 / 2)f(x - h) + (5 / 2)f(x + h) - 2f(x + 2h) + (1 / 2)f(x + 3h)) / h^5
    elseif 3 <= acc && acc <= 5
      # Accuracy 4
      return ((1 / 6)f(x - 4h) - (3 / 2)f(x - 3h) + (13 / 3)f(x - 2h) - (29 / 6)f(x - h) + (29 / 6)f(x + h) - (13 / 3)f(x + 2h) + (3 / 2)f(x + 3h) - (1 / 6)f(x + 4h)) / h^5
    elseif 5 <= acc
      # Accuracy 6
      return (-(13 / 288)f(x - 5h) + (19 / 36)f(x - 4h) - (87 / 32)f(x - 3h) + (13 / 2)f(x - 2h) - (323 / 48)f(x - h) + (323 / 48)f(x + h) - (13 / 2)f(x + 2h) + (87 / 32)f(x + 3h) - (19 / 36)f(x + 4h) + (13 / 288)f(x + 5h)) / h^5
    end
  elseif degree >= 6
    if acc <= 3
      # Accuracy 2
      return (f(x - 3h) - 6f(x - 2h) + 15f(x - h) - 20f(x) + 15f(x + h) - 6f(x + 2h) + f(x + 3h)) / h^6
    elseif 3 <= acc && acc <= 5
      # Accuracy 4
      return (-(1 / 4)f(x - 4h) + 3f(x - 3h) - 13f(x - 2h) + 29f(x - h) - (75 / 2)f(x) + 29f(x + h) - 13f(x + 2h) + 3f(x + 3h) - (1 / 4)f(x + 4h)) / h^6
    elseif 5 <= acc
      # Accuracy 6
      return ((13 / 240)f(x - 5h) - (19 / 24)f(x - 4h) + (87 / 16)f(x - 3h) - (39 / 2)f(x - 2h) + (323 / 8)f(x - h) - (1023 / 20)f(x) + (323 / 8)f(x - h) - (39 / 2)f(x - 2h) + (87 / 16)f(x - 3h) - (19 / 24)f(x - 4h) + (13 / 240)f(x - 5h)) / h^6
    end
  end
end

function __d_max_prec(degree::Integer, f::Function, x::Real, acc::Integer)::Real
  f_old = __d(degree, f, x, 0.1, acc)
  f_mid = __d(degree, f, x, 0.01, acc)
  f_new = __d(degree, f, x, 0.001, acc)
  δ_old = abs(f_old - f_mid)
  δ_new = abs(f_new - f_mid)
  h = 0.001
  while true
    h /= 10
    f_old = f_mid
    f_mid = f_new
    f_new = __d(degree, f, x, h, acc)
    δ_old = δ_new
    δ_new = abs(f_mid - f_new)
    if δ_new >= δ_old
      break
    end
  end
  return f_mid
end

function __d_h²(degree::Integer, f::Function, x::Real, acc::Integer, h::Float64=0.1, depth::Integer=4)::Real
  if depth > 9
    depth = 9
  end
  return __Dⱼₖ(degree, f, x, acc, h, 0, depth)
end

function __Dⱼₖ(degree::Integer, f::Function, x::Real, acc::Integer, h::Float64, i::Integer, k::Integer)::Real
  if k == 0
    return __d(degree, f, x, h / 2^i, acc)
  end
  return (4^k * __Dⱼₖ(degree, f, x, acc, h, i + 1, k - 1) - __Dⱼₖ(degree, f, x, acc, h, i, k - 1)) / (4^k - 1)
end

function __d_h²_max_prec(degree::Integer, f::Function, x::Real)::Real
  # Todo: this does not yield the highest precision
  function d_h²(degree::Integer, f::Function, x::Real, depth::Integer)::Real
    f_old = __d_h²(degree, f, x, 2, 0.1, depth)
    f_mid = __d_h²(degree, f, x, 2, 0.01, depth)
    f_new = __d_h²(degree, f, x, 2, 0.001, depth)
    δ_old = abs(f_old - f_mid)
    δ_new = abs(f_new - f_mid)
    h = 0.001
    while true
      h /= 10
      f_old = f_mid
      f_mid = f_new
      f_new = __d_h²(degree, f, x, 2, h, depth)
      δ_old = δ_new
      δ_new = abs(f_mid - f_new)
      if δ_new >= δ_old
        break
      end
    end
    return f_mid
  end

  # Todo: depth needs to be calculated differently
  f_old = d_h²(degree, f, x, 1)
  f_mid = d_h²(degree, f, x, 2)
  f_new = d_h²(degree, f, x, 3)
  δ_old = abs(f_old - f_mid)
  δ_new = abs(f_new - f_mid)
  if δ_new >= δ_old
    return f_mid
  end
  for depth in 3:7
    f_old = f_mid
    f_mid = f_new
    f_new = d_h²(degree, f, x, depth)
    δ_old = δ_new
    δ_new = abs(f_mid - f_new)
    if δ_new >= δ_old
      break
    end
  end
  return f_mid
end

# %%%%%%%%%%%%%%%%%%%% Partial derivative %%%%%%%%%%%%%%%%%%%%

"""
∂f/∂x
"""
function ∂(degree::Integer, f::Function, x⃗::Vector{<:Real}, i::Integer; max_prec=true, acc=8, ext=false, h=0.1, depth=4)::Real
  if ext
    # Todo: h² extrapolation
    if max_prec
      # Todo
    else
      # Todo
      return __∂_h²(degree, f, x⃗, acc, h, depth)
    end
  else
    if max_prec
      return __∂_max_prec(degree, f, x⃗, __d)
    else
      return __∂(degree, f, x⃗, i, h, acc)
    end
  end
end

function __∂(degree::Integer, f::Function, x⃗::Vector{<:Real}, i::Integer, h::Float64, acc::Integer)::Real
  # B/c accuracy is an integer, it's being clamped to avoid unnecessary exceptions
  h⃗ = zeros(length(x⃗))
  h⃗[i] = h
  if degree <= 1
    if acc <= 3
      # Accuracy 2
      return (-(1 / 2)f(x - h) + (1 / 2)f(x + h)) / h
    elseif 3 <= acc && acc <= 5
      # Accuracy 4
      return ((1 / 12)f(x - 2h) - (2 / 3)f(x - h) + (2 / 3)f(x + h) - (1 / 12)f(x + 2h)) / h
    elseif 5 <= acc && acc <= 7
      # Accuracy 6
      return (-(1 / 60)f(x - 3h) + (3 / 20)f(x - 2h) + (3 / 2)f(x - h) - (49 / 18)f(x) + (3 / 2)f(x + h) - (3 / 20)f(x + 2h) + (1 / 90)f(x + 3h)) / h
    elseif 7 <= acc
      # Accuracy 8
      return ((1 / 280)f(x - 4h) - (4 / 105)f(x - 3h) + (1 / 5)f(x - 2h) - (4 / 5)f(x - h) + (4 / 5)f(x + h) - (1 / 5)f(x + 2h) + (4 / 105)f(x + 3h) - (1 / 280)f(x + 4h)) / h
    end
  elseif degree == 2
    if acc <= 3
      # Accuracy 2
      return (f(x - h) - 2f(x) + f(x + h)) / h^2
    elseif 3 <= acc && acc <= 5
      # Accuracy 4
      return (-(1 / 12)f(x - 2h) + (4 / 3)f(x - h) - (5 / 2)f(x) + (4 / 3)f(x + h) - (1 / 12)f(x + 2h)) / h^2
    elseif 5 <= acc && acc <= 7
      # Accuracy 6
      return ((1 / 90)f(x - 3h) - (3 / 20)f(x - 2h) - (3 / 4)f(x - h) + (3 / 4)f(x + h) - (3 / 20)f(x + 2h) + (1 / 60)f(x + 3h)) / h^2
    elseif 7 <= acc
      # Accuracy 8
      return (-(1 / 560)f(x - 4h) + (8 / 315)f(x - 3h) - (1 / 5)f(x - 2h) + (8 / 5)f(x - h) - (205 / 72)f(x) + (8 / 5)f(x + h) - (1 / 5)f(x + 2h) + (8 / 315)f(x + 3h) - (1 / 560)f(x + 4h)) / h^2
    end
  elseif degree == 3
    if acc <= 3
      # Accuracy 2
      return (-(1 / 2)f(x - 2h) + f(x - h) - f(x + h) + (1 / 2)f(x + 2h)) / h^3
    elseif 3 <= acc && acc <= 5
      # Accuracy 4
      return ((1 / 8)f(x - 3h) - f(x - 2h) + (13 / 8)f(x - h) - (13 / 8)f(x + h) + f(x + 2h) - (1 / 8)f(x + 3h)) / h^3
    elseif 5 <= acc
      # Accuracy 6
      return (-(7 / 240)f(x - 4h) + (3 / 10)f(x - 3h) - (169 / 120)f(x - 2h) + (61 / 30)f(x - h) - (61 / 30)f(x + h) + (169 / 120)f(x + 2h) - (3 / 10)f(x + 3h) + (7 / 240)f(x + 4h)) / h^3
    end
  elseif degree == 4
    if acc <= 3
      # Accuracy 2
      return (f(x - 2h) - 4f(x - h) + 6f(x) - 4f(x + h) + f(x + 2h)) / h^4
    elseif 3 <= acc && acc <= 5
      # Accuracy 4
      return (-(1 / 6)f(x - 3h) + 2f(x - 2h) - (13 / 2)f(x - h) + (28 / 3)f(x) - (13 / 2)f(x + h) + 2f(x + 2h) - (1 / 6)f(x + 3h)) / h^4
    elseif 5 <= acc
      # Accuracy 6
      return ((7 / 240)f(x - 4h) - (2 / 5)f(x - 3h) + (169 / 60)f(x - 2h) - (122 / 15)f(x - h) + (91 / 8)f(x) - (122 / 15)f(x + h) + (169 / 60)f(x + 2h) - (2 / 5)f(x + 3h) + (7 / 240)f(x + 4h)) / h^4
    end
  elseif degree == 5
    if acc <= 3
      # Accuracy 2
      return (-(1 / 2)f(x - 3h) + 2f(x - 2h) - (5 / 2)f(x - h) + (5 / 2)f(x + h) - 2f(x + 2h) + (1 / 2)f(x + 3h)) / h^5
    elseif 3 <= acc && acc <= 5
      # Accuracy 4
      return ((1 / 6)f(x - 4h) - (3 / 2)f(x - 3h) + (13 / 3)f(x - 2h) - (29 / 6)f(x - h) + (29 / 6)f(x + h) - (13 / 3)f(x + 2h) + (3 / 2)f(x + 3h) - (1 / 6)f(x + 4h)) / h^5
    elseif 5 <= acc
      # Accuracy 6
      return (-(13 / 288)f(x - 5h) + (19 / 36)f(x - 4h) - (87 / 32)f(x - 3h) + (13 / 2)f(x - 2h) - (323 / 48)f(x - h) + (323 / 48)f(x + h) - (13 / 2)f(x + 2h) + (87 / 32)f(x + 3h) - (19 / 36)f(x + 4h) + (13 / 288)f(x + 5h)) / h^5
    end
  elseif degree >= 6
    if acc <= 3
      # Accuracy 2
      return (f(x - 3h) - 6f(x - 2h) + 15f(x - h) - 20f(x) + 15f(x + h) - 6f(x + 2h) + f(x + 3h)) / h^6
    elseif 3 <= acc && acc <= 5
      # Accuracy 4
      return (-(1 / 4)f(x - 4h) + 3f(x - 3h) - 13f(x - 2h) + 29f(x - h) - (75 / 2)f(x) + 29f(x + h) - 13f(x + 2h) + 3f(x + 3h) - (1 / 4)f(x + 4h)) / h^6
    elseif 5 <= acc
      # Accuracy 6
      return ((13 / 240)f(x - 5h) - (19 / 24)f(x - 4h) + (87 / 16)f(x - 3h) - (39 / 2)f(x - 2h) + (323 / 8)f(x - h) - (1023 / 20)f(x) + (323 / 8)f(x - h) - (39 / 2)f(x - 2h) + (87 / 16)f(x - 3h) - (19 / 24)f(x - 4h) + (13 / 240)f(x - 5h)) / h^6
    end
  end
end

function ∂(f::Function, x⃗::Vector{<:Real}, i::Integer; h::Float64=1e-4, max_prec::Bool=true)::Real
  if max_prec
    return __∂_max_prec(f, x⃗, i)
  end
  return __∂(f, x⃗, i, h)
end

function __∂(f::Function, x⃗::Vector{<:Real}, i::Integer, h::Float64)::Real
  h⃗ = zeros(length(x⃗))
  h⃗[i] = h
  return (-f(x⃗ + 2h⃗) + 8f(x⃗ + h⃗) - 8f(x⃗ - h⃗) + f(x⃗ - 2h⃗)) / 12h
end

function __∂_max_prec(f::Function, x⃗::Vector{<:Real}, i::Integer)::Real
  # Todo: rework this à la d()
  h = 10^-1
  f_old = __∂(f, x⃗, i, h)
  h /= 10
  f_new = __∂(f, x⃗, i, h)
  while true
    h /= 10
    f_old = f_new
    f_new = __∂(f, x⃗, i, h)
    f_old - f_new < 0 || break
  end
  return f_old
end

function ∂²(f::Function, x⃗::Vector{<:Real}, i::Integer; h::Float64=1e-4, max_prec::Bool=true)::Real
  if max_prec
    return __∂²_max_prec(f, x⃗, i)
  end
  return __∂²(f, x⃗, i, h)
end

function __∂²(f::Function, x⃗::Vector{<:Real}, i::Integer, h::Float64)::Real
  h⃗ = zeros(length(x⃗))
  h⃗[i] = h
  return (-f(x⃗ + 2h⃗) + 16f(x⃗ + h⃗) - 30f(x⃗) + 16f(x⃗ - h⃗) - f(x⃗ - 2h⃗)) / 12h^2
end

function __∂²_max_prec(f::Function, x⃗::Vector{<:Real}, i::Integer)::Real
  h = 10^-1
  f_old = __∂²(f, x⃗, i, h)
  h /= 10
  f_new = __∂²(f, x⃗, i, h)
  while true
    h /= 10
    f_old = f_new
    f_new = __∂²(f, x⃗, i, h)
    f_old - f_new < 0 || break
  end
  return f_old
end

function ∂³(f::Function, x⃗::Vector{<:Real}, i::Integer; h::Float64=1e-4, max_prec::Bool=true)::Real
  if max_prec
    return __∂³_max_prec(f, x⃗, i)
  end
  return __∂³(f, x⃗, i, h)
end

function __∂³(f::Function, x⃗::Vector{<:Real}, i::Integer, h::Float64)::Real
  h⃗ = zeros(length(x⃗))
  h⃗[i] = h
  return (f(x⃗ + 2h⃗) - 2f(x⃗ + h⃗) + 2f(x⃗ - h⃗) - f(x⃗ - 2h⃗)) / 2h^3
end

function __∂³_max_prec(f::Function, x⃗::Vector{<:Real}, i::Integer)::Real
  h = 10^-1
  f_old = __∂³(f, x⃗, i, h)
  h /= 10
  f_new = __∂³(f, x⃗, i, h)
  while true
    h /= 10
    f_old = f_new
    f_new = __∂³(f, x⃗, i, h)
    f_old - f_new < 0 || break
  end
  return f_old
end

function ∂⁴(f::Function, x⃗::Vector{<:Real}, i::Integer; h::Float64=1e-4, max_prec::Bool=true)::Real
  if max_prec
    return __∂⁴_max_prec(f, x⃗, i)
  end
  return __∂⁴(f, x⃗, i, h)
end

function __∂⁴(f::Function, x⃗::Vector{<:Real}, i::Integer, h::Float64)::Real
  h⃗ = zeros(length(x⃗))
  h⃗[i] = h
  return (f(x⃗ + 2h⃗) - 4f(x⃗ + h⃗) + 6f(x⃗) - 4f(x⃗ - h⃗) + f(x⃗ - 2h⃗)) / h^4
end

function __∂⁴_max_prec(f::Function, x⃗::Vector{<:Real}, i::Integer)::Real
  h = 10^-1
  f_old = __∂⁴(f, x⃗, i, h)
  h /= 10
  f_new = __∂⁴(f, x⃗, i, h)
  while true
    h /= 10
    f_old = f_new
    f_new = __∂⁴(f, x⃗, i, h)
    f_old - f_new < 0 || break
  end
  return f_old
end

# %%%%%%%%%%%%%%%%%%%% h² extrapolation for partials %%%%%%%%%%%%%%%%%%%%

function ∂_h²(f::Function, x⃗::Vector{Real}, i::Integer, h::Float64=1e-4)::Real
end

# %%%%%%%%%%%%%%%%%%%% Gradient %%%%%%%%%%%%%%%%%%%%

"""
Gradient
"""
function ∇(f::Function, x⃗::Vector{<:Real}; h::Float64=1e-4, max_prec::Bool=true)::Vector{<:Real}
  return [∂(f, x⃗, i; h=h, max_prec=max_prec) for i in 1:length(x⃗)]
end

function ∇²(f::Function, x⃗::Vector{<:Real}; h::Float64=1e-4, max_prec::Bool=true)::Vector{<:Real}
  return [∂²(f, x⃗, i; h=h, max_prec=max_prec) for i in 1:length(x⃗)]
end

function ∇³(f::Function, x⃗::Vector{<:Real}; h::Float64=1e-4, max_prec::Bool=true)::Vector{<:Real}
  return [∂³(f, x⃗, i; h=h, max_prec=max_prec) for i in 1:length(x⃗)]
end

function ∇⁴(f::Function, x⃗::Vector{<:Real}; h::Float64=1e-4, max_prec::Bool=true)::Vector{<:Real}
  return [∂⁴(f, x⃗, i; h=h, max_prec=max_prec) for i in 1:length(x⃗)]
end

# %%%%%%%%%%%%%%%%%%%% Jacobian %%%%%%%%%%%%%%%%%%%%

"""
Jacobian
"""
function D(f::Function, x⃗::Vector{<:Real}, h::Float64=1e-4)::Vector{<:Real}
end
