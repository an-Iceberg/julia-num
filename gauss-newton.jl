using LinearAlgebra

include("jacobian.jl")

"""
`f::(Real, Vector{<:Real})::Real` Function, which takes in a scalar `x` and an array of coefficients.

`λ₀` Initial coefficients

`x` x part of datapoints

`y` y part of datapoints

`h` Precision

`n` # of iterations to do

**Returns**:

Fitted coefficients

error

\\# of iterations
"""
function gauss_newton(f::Function, λ₀::Vector{<:Real}, x::Vector{<:Real}, y::Vector{<:Real}, h::Real=1e-3, n::Int=200)::Tuple{Vector{<:Real},Real,Int}
  g(λ::Vector{<:Real})::Vector{<:Real} = [y - f(x, λ) for (x, y) in zip(x, y)]
  Ẽ(λ::Vector{<:Real})::Real = norm(g(λ))^2
  increment = h + 1
  iter_count = 0
  λ = copy(λ₀)

  # Todo: implement
  while iter_count < n && increment > h
    Q, R = qr(D1_4(g, λ, h))
    break
  end

  return (λ, Ẽ(λ), iter_count)
end

x = [0, 1, 2, 3, 4]
y = [3, 1, 0.5, 0.2, 0.05]

function f(x::Real, coefs::Vector{<:Real})::Real
  a, b = coefs
  return a * ℯ^(b * x)
end

h = 1e-7
λ = [2.0, 2.0]

λ, Ẽ, n = gauss_newton(f, λ, x, y, h)
a = round(λ[1]; digits=2)
b = round(λ[2]; digits=2)
println("f(x) = $a⋅ℯ^($b⋅x)")
println("Ẽ = $Ẽ")
