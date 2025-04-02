function f(x::Real, coefs::Vector{<:Real})::Real
  a, b = coefs
  return a * x + b
end

function gauss_newton(f::Function, λ₀::Vector{<:Real}, ϵ::Real, n::Integer)::Tuple{Vector{<:Real},Integer}
end
