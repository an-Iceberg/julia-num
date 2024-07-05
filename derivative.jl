# %%%%%%%%%%%%%%%%%%%% Midpoint derivative %%%%%%%%%%%%%%%%%%%%

function __d_midpoint(f::Function, x::Real, h::Float64=1e-4)::Real
  return (f(x + h) - f(x - h)) / 2h
end

function __d²_midpoint(f::Function, x::Real, h::Float64=1e-4)::Real
  return (f(x + h) - 2f(x) + f(x - h)) / h^2
end

function __d³_midpoint(f::Function, x::Real, h::Float64=1e-4)::Real
  return (-f(x - 2h) + 2f(x - h) - 2f(x + h) + f(x + 2h)) / 2h^3
end

function __d⁴_midpoint(f::Function, x::Real, h::Float64=1e-4)::Real
  return (f(x - 2h) - 4f(x - h) + 6f(x) - 4f(x + h) + f(x + 2h)) / h^4
end

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

# Todo: make method (5 point, h²) choosable via string
# Todo: make options (method, max_prec, h, depth) an optional tuple with default values
"""
df/dx

Calculates the first derivative of a function at a point numerically.

---

# Parameters
`f (Real) -> Real` The function in question.

`x` The point at which derivation takes place.

## Optional
`h` Precision.

`max_prec` Whether to find the value with the most precision. Enabling this will
take a little longer to compute the result, but it will have the highest precision.

# Returns
The derivative of `f` at `x`.
"""
function d(f::Function, x::Real; h::Float64=1e-4, max_prec::Bool=true)::Real
  if max_prec
    return __d_max_prec(f, x)
  end
  return __d(f, x, h)
end

function __d(f::Function, x::Real, h::Float64)::Real
  return (-f(x + 2h) + 8f(x + h) - 8f(x - h) + f(x - 2h)) / 12h
end

function __d_max_prec(f::Function, x::Real)::Real
  h = 1e-1
  f_old = 0
  f_mid = 0
  f_new = __d(f, x, h)
  δ_old = 0
  δ_new = Inf64
  while true
    h /= 5
    f_old = f_mid
    f_mid = f_new
    f_new = __d(f, x, h)
    δ_old = δ_new
    δ_new = abs(f_mid - f_new)
    if δ_new >= δ_old
      break
    end
  end
  return f_old
end

function d²(f::Function, x::Real; h::Float64=1e-4, max_prec::Bool=true)::Real
  if max_prec
    return __d²_max_prec(f, x)
  end
  return __d²(f, x, h)
end

function __d²(f::Function, x::Real, h::Float64)::Real
  return (-f(x + 2h) + 16f(x + h) - 30f(x) + 16f(x - h) - f(x - 2h)) / 12h^2
end

function __d²_max_prec(f::Function, x::Real)::Real
  h = 1e-1
  f_old = 0
  f_mid = 0
  f_new = __d²(f, x, h)
  δ_old = 0
  δ_new = Inf64
  while true
    h /= 5
    f_old = f_mid
    f_mid = f_new
    f_new = __d²(f, x, h)
    δ_old = δ_new
    δ_new = abs(f_mid - f_new)
    if δ_new >= δ_old
      break
    end
  end
  return f_old
end

function d³(f::Function, x::Real; h::Float64=1e-4, max_prec::Bool=true)::Real
  if max_prec
    return __d³_max_prec(f, x)
  end
  return __d³(f, x, h)
end

function __d³(f::Function, x::Real, h::Float64)::Real
  return (f(x + 2h) - 2f(x + h) + 2f(x - h) - f(x - 2h)) / 2h^3
end

function __d³_max_prec(f::Function, x::Real)::Real
  h = 1e-1
  f_old = 0
  f_mid = 0
  f_new = __d³(f, x, h)
  δ_old = 0
  δ_new = Inf64
  while true
    h /= 5
    f_old = f_mid
    f_mid = f_new
    f_new = __d³(f, x, h)
    δ_old = δ_new
    δ_new = abs(f_mid - f_new)
    if δ_new >= δ_old
      break
    end
  end
  return f_old
end

function d⁴(f::Function, x::Real; h::Float64=1e-4, max_prec::Bool=true)::Real
  if max_prec
    return __d⁴_max_prec(f, x)
  end
  return __d⁴(f, x, h)
end

function __d⁴(f::Function, x::Real, h::Float64)::Real
  return (f(x + 2h) - 4f(x + h) + 6f(x) - 4f(x - h) + f(x - 2h)) / h^4
end

function __d⁴_max_prec(f::Function, x::Real)::Real
  h = 1e-1
  f_old = 0
  f_mid = 0
  f_new = __d⁴(f, x, h)
  δ_old = 0
  δ_new = Inf64
  while true
    h /= 5
    f_old = f_mid
    f_mid = f_new
    f_new = __d⁴(f, x, h)
    δ_old = δ_new
    δ_new = abs(f_mid - f_new)
    if δ_new >= δ_old
      break
    end
  end
  return f_old
end

# %%%%%%%%%%%%%%%%%%%% h² extrapolation %%%%%%%%%%%%%%%%%%%%

function __d_h²(f::Function, x::Real, h::Float64=0.1, depth::Integer=4)::Real
  if depth > 9
    depth = 9
  end
  return __Dⱼₖ(f, x, h, 0, depth)
end

function __Dⱼₖ(f::Function, x::Real, h::Float64, i::Integer, k::Integer)::Real
  if k == 0
    return __d_midpoint(f, x, h / 2^i)
  end
  return (4^k * __Dⱼₖ(f, x, h, i + 1, k - 1) - __Dⱼₖ(f, x, h, i, k - 1)) / (4^k - 1)
end

function __d_h²_max_prec(f::Function, x::Real)::Real

  function d_h²(f::Function, x::Real, depth::Integer)::Real
    h = 1e-1
    f_old = 0
    f_mid = 0
    f_new = __d_h²(f, x, h)
    δ_old = 0
    δ_new = Inf64
    while true
      h /= 5
      f_old = f_mid
      f_mid = f_new
      f_new = __d_h²(f, x, h)
      δ_old = δ_new
      δ_new = abs(f_mid - f_new)
      if δ_new >= δ_old
        break
      end
    end
    return f_old
  end

  f_old = 0
  f_mid = 0
  f_new = d_h²(f, x, 1)
  δ_old = 0
  δ_new = Inf64
  for depth in 2:7
    f_old = f_mid
    f_mid = f_new
    f_new = d_h²(f, x, depth)
    δ_old = δ_new
    δ_new = abs(f_mid - f_new)
    if δ_new >= δ_old
      break
    end
  end
  return f_old
end

function __d²_h²(f::Function, x::Real, h::Float64=0.1, depth::Integer=4)::Real
  if depth > 9
    depth = 9
  end
  return __D²ⱼₖ(f, x, h, 0, depth)
end

function __D²ⱼₖ(f::Function, x::Real, h::Float64, i::Integer, k::Integer)::Real
  if k == 0
    return __d²_midpoint(f, x, h / 2^i)
  end
  return (4^k * __D²ⱼₖ(f, x, h, i + 1, k - 1) - __D²ⱼₖ(f, x, h, i, k - 1)) / (4^k - 1)
end

function __d²_h²_max_prec(f::Function, x::Real)::Real

  function d²_h²(f::Function, x::Real, depth::Integer)::Real
    h = 1e-1
    f_old = 0
    f_mid = 0
    f_new = __d²_h²(f, x, h)
    δ_old = 0
    δ_new = Inf64
    while true
      h /= 5
      f_old = f_mid
      f_mid = f_new
      f_new = __d²_h²(f, x, h)
      δ_old = δ_new
      δ_new = abs(f_mid - f_new)
      if δ_new >= δ_old
        break
      end
    end
    return f_old
  end

  f_old = 0
  f_mid = 0
  f_new = d²_h²(f, x, 1)
  δ_old = 0
  δ_new = Inf64
  for depth in 2:7
    f_old = f_mid
    f_mid = f_new
    f_new = d²_h²(f, x, depth)
    δ_old = δ_new
    δ_new = abs(f_mid - f_new)
    if δ_new >= δ_old
      break
    end
  end
  return f_old
end

function __d³_h²(f::Function, x::Real, h::Float64=0.1, depth::Integer=4)::Real
  if depth > 9
    depth = 9
  end
  return __D³ⱼₖ(f, x, h, 0, depth)
end

function __D³ⱼₖ(f::Function, x::Real, h::Float64, i::Integer, k::Integer)::Real
  if k == 0
    return __d³_midpoint(f, x, h / 2^i)
  end
  return (4^k * __D³ⱼₖ(f, x, h, i + 1, k - 1) - __D³ⱼₖ(f, x, h, i, k - 1)) / (4^k - 1)
end

function __d³_h²_max_prec(f::Function, x::Real)::Real

  function d³_h²(f::Function, x::Real, depth::Integer)::Real
    h = 1e-1
    f_old = 0
    f_mid = 0
    f_new = __d³_h²(f, x, h)
    δ_old = 0
    δ_new = Inf64
    while true
      h /= 5
      f_old = f_mid
      f_mid = f_new
      f_new = __d³_h²(f, x, h)
      δ_old = δ_new
      δ_new = abs(f_mid - f_new)
      if δ_new >= δ_old
        break
      end
    end
    return f_old
  end

  f_old = 0
  f_mid = 0
  f_new = d³_h²(f, x, 1)
  δ_old = 0
  δ_new = Inf64
  for depth in 2:7
    f_old = f_mid
    f_mid = f_new
    f_new = d³_h²(f, x, depth)
    δ_old = δ_new
    δ_new = abs(f_mid - f_new)
    if δ_new >= δ_old
      break
    end
  end
  return f_old
end

function __d⁴_h²(f::Function, x::Real, h::Float64=0.1, depth::Integer=4)::Real
  if depth > 9
    depth = 9
  end
  return __D⁴ⱼₖ(f, x, h, 0, depth)
end

function __D⁴ⱼₖ(f::Function, x::Real, h::Float64, i::Integer, k::Integer)::Real
  if k == 0
    return __d⁴_midpoint(f, x, h / 2^i)
  end
  return (4^k * __D⁴ⱼₖ(f, x, h, i + 1, k - 1) - __D⁴ⱼₖ(f, x, h, i, k - 1)) / (4^k - 1)
end

function __d⁴_h²_max_prec(f::Function, x::Real)::Real

  function d⁴_h²(f::Function, x::Real, depth::Integer)::Real
    h = 1e-1
    f_old = 0
    f_mid = 0
    f_new = __d⁴_h²(f, x, h)
    δ_old = 0
    δ_new = Inf64
    while true
      h /= 5
      f_old = f_mid
      f_mid = f_new
      f_new = __d⁴_h²(f, x, h)
      δ_old = δ_new
      δ_new = abs(f_mid - f_new)
      if δ_new >= δ_old
        break
      end
    end
    return f_old
  end

  f_old = 0
  f_mid = 0
  f_new = d⁴_h²(f, x, 1)
  δ_old = 0
  δ_new = Inf64
  for depth in 2:7
    f_old = f_mid
    f_mid = f_new
    f_new = d⁴_h²(f, x, depth)
    δ_old = δ_new
    δ_new = abs(f_mid - f_new)
    if δ_new >= δ_old
      break
    end
  end
  return f_old
end

# %%%%%%%%%%%%%%%%%%%% Partial derivative %%%%%%%%%%%%%%%%%%%%

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
