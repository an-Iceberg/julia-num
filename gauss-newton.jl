function f(x::Real, coefs::Vector{<:Real})::Real
  a, b = coefs
  return a * x + b
end
