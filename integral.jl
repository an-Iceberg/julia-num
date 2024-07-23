# Todo: restructure this
function ∫(f::Function, a::Real, b::Real; max_prec::Bool=true, ext::Bool=false, h::Float64=0.1, depth::Integer=4)::Real
  if ext
    if max_prec
    else
      return __∫_romberg(f, a, b, h, depth)
    end
  else
    if max_prec
    else
      return __∫(f, a, b, h)
    end
  end
end

# %%%%%%%%%%%%%%%%%%%% Simpson's ⅓ rule %%%%%%%%%%%%%%%%%%%%

function __∫_max_prec(f::Function, a::Real, b::Real)::Real
end

function __∫(f::Function, a::Real, b::Real, h::Float64)::Real
  # Integral using Simpson's rule.
  n = Int(ceil((b - a) / h))
  Σ1 = 0
  Σ2 = 0
  x(i) = a + i * h
  for i in 1:(n-1)
    Σ1 += f(x(i))
  end
  for i in 1:n
    Σ2 += f((x(i - 1) + x(i)) / 2)
  end
  return (h / 3) * (0.5f(a) + Σ1 + 2Σ2 + 0.5f(b))
end

# %%%%%%%%%%%%%%%%%%%% Boole's rule %%%%%%%%%%%%%%%%%%%%
# https://en.wikipedia.org/wiki/Boole%27s_rule

# %%%%%%%%%%%%%%%%%%%% Romberg extrapolation %%%%%%%%%%%%%%%%%%%%

function __∫_romberg_max_prec(f::Function, a::Real, b::Real)::Real
end

function __∫_romberg(f::Function, a::Real, b::Real, h::Float64=0.1, depth::Integer=4)::Real
  if depth > 9
    depth = 9
  end
  return __Tⱼₖ(f, a, b, h, 0, depth)
end

function __Tⱼₖ(f::Function, a::Real, b::Real, h::Float64, j::Integer, k::Integer)::Real
  if k == 0
    return __T(f, a, b, h / 2^j)
  end
  return (4^k * __Tⱼₖ(f, a, b, h, j + 1, k - 1) - __Tⱼₖ(f, a, b, h, j, k - 1)) / (4^k - 1)
end

function __T(f::Function, a::Real, b::Real, h::Float64)::Real
  n = Int(ceil((b - a) / h))
  Σ = 0
  x(i) = a + i * h
  for i in 1:n
    Σ += f(x(i))
  end
  return h * (((f(a) + f(b)) / 2) + Σ)
end
