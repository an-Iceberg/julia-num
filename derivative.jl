# https://en.wikipedia.org/wiki/Finite_difference_coefficient

# %%%%%%%%%%%%%%%%%%%% Derivative %%%%%%%%%%%%%%%%%%%%

"""
df/dx

`f::(Real) -> Real`
"""
function d(degree::Integer, f::Function, x::Real; max_prec::Bool=true, acc::Integer=8, ext::Bool=false, h::Real=0.1, depth::Integer=4)::Real
  if ext
    if max_prec
      return __d_h²_max_prec(degree, f, x)
    else
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

function __d(degree::Integer, f::Function, x::Real, h::Real, acc::Integer)::Real
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
      return (-(1 / 60)f(x - 3h) + 0.15f(x - 2h) - 0.75f(x - h) + 0.75f(x + h) - 0.15f(x + 2h) + (1 / 60)f(x + 3h)) / h
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
  # ToDo: why does this not work?!?!?!
  f_1 = __d(degree, f, x, 1.0, acc)
  f_old = __d(degree, f, x, 0.1, acc)
  f_new = __d(degree, f, x, 0.01, acc)
  δ_old = abs(f_1 - f_old)
  println("    δ = ", δ_old)
  δ_new = abs(f_1 - f_new)
  println("    δ = ", δ_new)
  # In case f_1 is the most precise
  # if δ_δ > 0
  #   return f_1
  # end
  h = 0.01
  while δ_new > δ_old
    h /= 10
    f_old = f_new
    f_new = __d(degree, f, x, h, acc)
    δ_old = δ_new
    δ_new = abs(f_1 - f_new)
    println("    δ = ", δ_new)
  end
  return f_old
end

function __d_h²(degree::Integer, f::Function, x::Real, acc::Integer, h::Real=0.1, depth::Integer=4)::Real
  if depth > 9
    depth = 9
  end
  return __Dⱼₖ(degree, f, x, acc, h, 0, depth)
end

function __Dⱼₖ(degree::Integer, f::Function, x::Real, acc::Integer, h::Real, j::Integer, k::Integer)::Real
  if k == 0
    return __d(degree, f, x, h / 2^j, acc)
  end
  return (4^k * __Dⱼₖ(degree, f, x, acc, h, j + 1, k - 1) - __Dⱼₖ(degree, f, x, acc, h, j, k - 1)) / (4^k - 1)
end

# Todo: this does not yield the highest precision
# Idea: use gradient descent (x axis = h, y axis = depth)
function __d_h²_max_prec(degree::Integer, f::Function, x::Real)::Real

  function d_h²(degree::Integer, f::Function, x::Real, depth::Integer)::Real
    f_old = __d_h²(degree, f, x, 2, 1.0, depth)
    f_mid = __d_h²(degree, f, x, 2, 0.1, depth)
    f_new = __d_h²(degree, f, x, 2, 0.01, depth)
    δ_old = abs(f_old - f_mid)
    δ_new = abs(f_new - f_mid)
    h = 0.01
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
  f_old = d_h²(degree, f, x, 0)
  f_mid = d_h²(degree, f, x, 1)
  f_new = d_h²(degree, f, x, 2)
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

`f::(Vector{Real}) -> Real`
"""
function ∂(degree::Integer, f::Function, x⃗::Vector{<:Real}, i::Integer; max_prec::Bool=true, acc::Integer=8, ext::Bool=false, h::Real=0.1, depth::Integer=4)::Real
  if ext
    if max_prec
      return __∂_h²_max_prec(degree, f, x⃗, i)
    else
      return __∂_h²(degree, f, x⃗, i, acc, h, depth)
    end
  else
    if max_prec
      return __∂_max_prec(degree, f, x⃗, i, acc)
    else
      return __∂(degree, f, x⃗, i, h, acc)
    end
  end
end

function __∂(degree::Integer, f::Function, x⃗::Vector{<:Real}, i::Integer, h::Real, acc::Integer)::Real
  # B/c accuracy is an integer, it's being clamped to avoid unnecessary exceptions
  h⃗ = zeros(length(x⃗))
  h⃗[i] = h
  if degree <= 1
    if acc <= 3
      # Accuracy 2
      return (-(1 / 2)f(x⃗ - h⃗) + (1 / 2)f(x⃗ + h⃗)) / h
    elseif 3 <= acc && acc <= 5
      # Accuracy 4
      return ((1 / 12)f(x⃗ - 2h⃗) - (2 / 3)f(x⃗ - h⃗) + (2 / 3)f(x⃗ + h⃗) - (1 / 12)f(x⃗ + 2h⃗)) / h
    elseif 5 <= acc && acc <= 7
      # Accuracy 6
      return (-(1 / 60)f(x⃗ - 3h⃗) + (3 / 20)f(x⃗ - 2h⃗) + (3 / 2)f(x⃗ - h⃗) - (49 / 18)f(x⃗) + (3 / 2)f(x⃗ + h⃗) - (3 / 20)f(x⃗ + 2h⃗) + (1 / 60)f(x⃗ + 3h⃗)) / h
    elseif 7 <= acc
      # Accuracy 8
      return ((1 / 280)f(x⃗ - 4h⃗) - (4 / 105)f(x⃗ - 3h⃗) + (1 / 5)f(x⃗ - 2h⃗) - (4 / 5)f(x⃗ - h⃗) + (4 / 5)f(x⃗ + h⃗) - (1 / 5)f(x⃗ + 2h⃗) + (4 / 105)f(x⃗ + 3h⃗) - (1 / 280)f(x⃗ + 4h⃗)) / h
    end
  elseif degree == 2
    if acc <= 3
      # Accuracy 2
      return (f(x⃗ - h⃗) - 2f(x⃗) + f(x⃗ + h⃗)) / h^2
    elseif 3 <= acc && acc <= 5
      # Accuracy 4
      return (-(1 / 12)f(x⃗ - 2h⃗) + (4 / 3)f(x⃗ - h⃗) - (5 / 2)f(x⃗) + (4 / 3)f(x⃗ + h⃗) - (1 / 12)f(x⃗ + 2h⃗)) / h^2
    elseif 5 <= acc && acc <= 7
      # Accuracy 6
      return ((1 / 90)f(x⃗ - 3h⃗) - (3 / 20)f(x⃗ - 2h⃗) + (3 / 2)f(x⃗ - h⃗) - (49 / 18)f(x⃗) + (3 / 2)f(x⃗ + h⃗) - (3 / 20)f(x⃗ + 2h⃗) + (1 / 90)f(x⃗ + 3h⃗)) / h^2
    elseif 7 <= acc
      # Accuracy 8
      return (-(1 / 560)f(x⃗ - 4h⃗) + (8 / 315)f(x⃗ - 3h⃗) - (1 / 5)f(x⃗ - 2h⃗) + (8 / 5)f(x⃗ - h⃗) - (205 / 72)f(x⃗) + (8 / 5)f(x⃗ + h⃗) - (1 / 5)f(x⃗ + 2h⃗) + (8 / 315)f(x⃗ + 3h⃗) - (1 / 560)f(x⃗ + 4h⃗)) / h^2
    end
  elseif degree == 3
    if acc <= 3
      # Accuracy 2
      return (-(1 / 2)f(x⃗ - 2h⃗) + f(x⃗ - h⃗) - f(x⃗ + h⃗) + (1 / 2)f(x⃗ + 2h⃗)) / h^3
    elseif 3 <= acc && acc <= 5
      # Accuracy 4
      return ((1 / 8)f(x⃗ - 3h⃗) - f(x⃗ - 2h⃗) + (13 / 8)f(x⃗ - h⃗) - (13 / 8)f(x⃗ + h⃗) + f(x⃗ + 2h⃗) - (1 / 8)f(x⃗ + 3h⃗)) / h^3
    elseif 5 <= acc
      # Accuracy 6
      return (-(7 / 240)f(x⃗ - 4h⃗) + (3 / 10)f(x⃗ - 3h⃗) - (169 / 120)f(x⃗ - 2h⃗) + (61 / 30)f(x⃗ - h⃗) - (61 / 30)f(x⃗ + h⃗) + (169 / 120)f(x⃗ + 2h⃗) - (3 / 10)f(x⃗ + 3h⃗) + (7 / 240)f(x⃗ + 4h⃗)) / h^3
    end
  elseif degree == 4
    if acc <= 3
      # Accuracy 2
      return (f(x⃗ - 2h⃗) - 4f(x⃗ - h⃗) + 6f(x⃗) - 4f(x⃗ + h⃗) + f(x⃗ + 2h⃗)) / h^4
    elseif 3 <= acc && acc <= 5
      # Accuracy 4
      return (-(1 / 6)f(x⃗ - 3h⃗) + 2f(x⃗ - 2h⃗) - (13 / 2)f(x⃗ - h⃗) + (28 / 3)f(x⃗) - (13 / 2)f(x⃗ + h⃗) + 2f(x⃗ + 2h⃗) - (1 / 6)f(x⃗ + 3h⃗)) / h^4
    elseif 5 <= acc
      # Accuracy 6
      return ((7 / 240)f(x⃗ - 4h⃗) - (2 / 5)f(x⃗ - 3h⃗) + (169 / 60)f(x⃗ - 2h⃗) - (122 / 15)f(x⃗ - h⃗) + (91 / 8)f(x⃗) - (122 / 15)f(x⃗ + h⃗) + (169 / 60)f(x⃗ + 2h⃗) - (2 / 5)f(x⃗ + 3h⃗) + (7 / 240)f(x⃗ + 4h⃗)) / h^4
    end
  elseif degree == 5
    if acc <= 3
      # Accuracy 2
      return (-(1 / 2)f(x⃗ - 3h⃗) + 2f(x⃗ - 2h⃗) - (5 / 2)f(x⃗ - h⃗) + (5 / 2)f(x⃗ + h⃗) - 2f(x⃗ + 2h⃗) + (1 / 2)f(x⃗ + 3h⃗)) / h^5
    elseif 3 <= acc && acc <= 5
      # Accuracy 4
      return ((1 / 6)f(x⃗ - 4h⃗) - (3 / 2)f(x⃗ - 3h⃗) + (13 / 3)f(x⃗ - 2h⃗) - (29 / 6)f(x⃗ - h⃗) + (29 / 6)f(x⃗ + h⃗) - (13 / 3)f(x⃗ + 2h⃗) + (3 / 2)f(x⃗ + 3h⃗) - (1 / 6)f(x⃗ + 4h⃗)) / h^5
    elseif 5 <= acc
      # Accuracy 6
      return (-(13 / 288)f(x⃗ - 5h⃗) + (19 / 36)f(x⃗ - 4h⃗) - (87 / 32)f(x⃗ - 3h⃗) + (13 / 2)f(x⃗ - 2h⃗) - (323 / 48)f(x⃗ - h⃗) + (323 / 48)f(x⃗ + h⃗) - (13 / 2)f(x⃗ + 2h⃗) + (87 / 32)f(x⃗ + 3h⃗) - (19 / 36)f(x⃗ + 4h⃗) + (13 / 288)f(x⃗ + 5h⃗)) / h^5
    end
  elseif degree >= 6
    if acc <= 3
      # Accuracy 2
      return (f(x⃗ - 3h⃗) - 6f(x⃗ - 2h⃗) + 15f(x⃗ - h⃗) - 20f(x⃗) + 15f(x⃗ + h⃗) - 6f(x⃗ + 2h⃗) + f(x⃗ + 3h⃗)) / h^6
    elseif 3 <= acc && acc <= 5
      # Accuracy 4
      return (-(1 / 4)f(x⃗ - 4h⃗) + 3f(x⃗ - 3h⃗) - 13f(x⃗ - 2h⃗) + 29f(x⃗ - h⃗) - (75 / 2)f(x⃗) + 29f(x⃗ + h⃗) - 13f(x⃗ + 2h⃗) + 3f(x⃗ + 3h⃗) - (1 / 4)f(x⃗ + 4h⃗)) / h^6
    elseif 5 <= acc
      # Accuracy 6
      return ((13 / 240)f(x⃗ - 5h⃗) - (19 / 24)f(x⃗ - 4h⃗) + (87 / 16)f(x⃗ - 3h⃗) - (39 / 2)f(x⃗ - 2h⃗) + (323 / 8)f(x⃗ - h⃗) - (1023 / 20)f(x⃗) + (323 / 8)f(x⃗ - h⃗) - (39 / 2)f(x⃗ - 2h⃗) + (87 / 16)f(x⃗ - 3h⃗) - (19 / 24)f(x⃗ - 4h⃗) + (13 / 240)f(x⃗ - 5h⃗)) / h^6
    end
  end
end

function __∂_max_prec(degree::Integer, f::Function, x⃗::Vector{<:Real}, i::Integer, acc::Integer)::Real
  f_old = __∂(degree, f, x⃗, i, 0.1, acc)
  f_mid = __∂(degree, f, x⃗, i, 0.01, acc)
  f_new = __∂(degree, f, x⃗, i, 0.001, acc)
  δ_old = abs(f_old - f_mid)
  δ_new = abs(f_new - f_mid)
  h = 0.001
  while true
    h /= 10
    f_old = f_mid
    f_mid = f_new
    f_new = __∂(degree, f, x⃗, i, h, acc)
    δ_old = δ_new
    δ_new = abs(f_mid - f_new)
    if δ_new >= δ_old
      break
    end
  end
  return f_mid
end

function __∂_h²(degree::Integer, f::Function, x⃗::Vector{<:Real}, i::Integer, acc::Integer, h::Real=0.1, depth::Integer=4)::Real
  if depth > 9
    depth = 9
  end
  return __∂ⱼₖ(degree, f, x⃗, i, acc, h, 0, depth)
end

function __∂ⱼₖ(degree::Integer, f::Function, x⃗::Vector{<:Real}, i::Integer, acc::Integer, h::Real, j::Integer, k::Integer)::Real
  if k == 0
    return __∂(degree, f, x⃗, i, h / 2^j, acc)
  end
  return (4^k * __∂ⱼₖ(degree, f, x⃗, i, acc, h, j + 1, k - 1) - __∂ⱼₖ(degree, f, x⃗, i, acc, h, j, k - 1)) / (4^k - 1)
end

# Todo: this does not yield the highest precision
function __∂_h²_max_prec(degree::Integer, f::Function, x⃗::Vector{<:Real}, i::Integer)::Real

  function ∂_h²(degree::Integer, f::Function, x⃗::Vector{<:Real}, i::Integer, depth::Integer)::Real
    f_old = __∂_h²(degree, f, x⃗, i, 2, 0.1, depth)
    f_mid = __∂_h²(degree, f, x⃗, i, 2, 0.01, depth)
    f_new = __∂_h²(degree, f, x⃗, i, 2, 0.001, depth)
    δ_old = abs(f_old - f_mid)
    δ_new = abs(f_new - f_mid)
    h = 0.001
    while true
      h /= 10
      f_old = f_mid
      f_mid = f_new
      f_new = __∂_h²(degree, f, x⃗, i, 2, h, depth)
      δ_old = δ_new
      δ_new = abs(f_mid - f_new)
      if δ_new >= δ_old
        break
      end
    end
    return f_mid
  end

  # Todo: depth needs to be calculated differently
  f_old = ∂_h²(degree, f, x⃗, i, 1)
  f_mid = ∂_h²(degree, f, x⃗, i, 2)
  f_new = ∂_h²(degree, f, x⃗, i, 3)
  δ_old = abs(f_old - f_mid)
  δ_new = abs(f_new - f_mid)
  if δ_new >= δ_old
    return f_mid
  end
  for depth in 4:7
    f_old = f_mid
    f_mid = f_new
    f_new = ∂_h²(degree, f, x⃗, i, depth)
    δ_old = δ_new
    δ_new = abs(f_mid - f_new)
    if δ_new >= δ_old
      break
    end
  end
  return f_mid
end

# %%%%%%%%%%%%%%%%%%%% Gradient %%%%%%%%%%%%%%%%%%%%

"""
Gradient

`f::(Vector{Real}) -> Vector{Real}`
"""
function ∇(degree::Integer, f::Function, x⃗::Vector{<:Real}; max_prec::Bool=true, acc::Integer=8, ext::Bool=false, h::Real=0.1, depth::Integer=4)::Vector{<:Real}
  return [∂(degree, f, x⃗, i; max_prec=max_prec, acc=acc, ext=ext, h=h, depth=depth) for i in 1:length(x⃗)]
end

# %%%%%%%%%%%%%%%%%%%% Jacobian %%%%%%%%%%%%%%%%%%%%

# Todo: ∂ for f::(Vector{<:Real}) -> Vector{<:Real}

function ∂⃗(degree::Integer, f::Function, x⃗::Vector{<:Real}, i::Integer; max_prec::Bool=true, acc::Integer=8, ext::Bool=false, h::Real, depth::Integer=4)::Vector{<:Real}
end

function __∂⃗(degree::Integer, f::Function, x⃗::Vector{<:Real}, i::Integer, h::Real, acc::Integer)::Vector{<:Real}
  h⃗ = zeros(length(x⃗))
  h⃗[i] = h
  if degree <= 1
    if acc <= 3
      # Accuracy 2
      # ToDo: confirm that this works correctly
      return [n / h for n in -(1 / 2)f(x⃗ - h⃗) + (1 / 2)f(x⃗ + h⃗)]
    end
  elseif degree == 2
  elseif degree == 3
  elseif degree == 4
  elseif degree == 5
  elseif degree <= 6
  end
end

"""
Jacobian

`f::(Vector{Real}) -> Vector{Real}`
"""
function D(degree::Integer, f::Function, x⃗::Vector{<:Real}; max_prec::Bool=true, acc::Integer=8, ext::Bool=false, h::Real=0.1, depth::Integer=4)::Vector{<:Real}
  sol = []
  for i in 1:length(x⃗)
  end
end
