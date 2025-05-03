using ForwardDiff, LinearAlgebra

"""
`f::(Real, Vector{<:Real})::Real` Function, which takes in a scalar `x` and an array of coefficients.

`λ₀` Initial coefficients

`x` x part of datapoints

`y` y part of datapoints

`h` Precision

`n` # of iterations to do

**Returns**:

(Fitted coefficients, mean square error, \\# of iterations)
"""
function fit(
  f::Function,
  λ0::Vector{<:Number},
  x::Vector{<:Number},
  y::Vector{<:Number},
  h::Real=1e-3,
  n::Int=100,
)::Tuple{Vector{<:Number}, Real, Int}
  D = ForwardDiff.jacobian
  g(λ::Vector{<:Real})::Vector{<:Real} = [y - f(x, λ) for (x, y) in zip(x, y)]
  Ẽ(λ::Vector{<:Real})::Real = norm(g(λ))^2

  inc = h + 1 # Increment
  iter = 0 # Iteration count
  λ = copy(λ0)

  while iter < n && inc > h
    Q, R = qr(D(g, λ))
    δ = R \ (.-transpose(Matrix(Q)) * g(λ))
    λ += δ
    inc = norm(δ)
    iter += 1
  end

  return (λ, Ẽ(λ), iter)
end

"""
`f::(Real, Vector{<:Real})::Real` Function, which takes in a scalar `x` and an array of coefficients.

`λ₀` Initial coefficients

`x` x part of datapoints

`y` y part of datapoints

`h` Precision

`n` # of iterations to do

**Returns**:

(Fitted coefficients, mean square error, \\# of iterations)
"""
function fit_damped(
  f::Function,
  λ₀::Vector{<:Number},
  x::Vector{<:Number},
  y::Vector{<:Number},
  h::Real=1e-3,
  n::Int=100,
)::Tuple{Vector{<:Number}, Real, Int}
  D = ForwardDiff.jacobian
  g(λ::Vector{<:Real})::Vector{<:Real} = [y - f(x, λ) for (x, y) in zip(x, y)]
  Ẽ(λ::Vector{<:Real})::Real = norm(g(λ))^2

  inc = h + 1 # Increment
  iter = 0 # Iteration count
  λ = copy(λ₀)

  while iter < n && inc > h
    Q, R = qr(D(g, λ))
    δ = R \ (.-transpose(Matrix(Q)) * g(λ))

    # Calculating how strong the damping should be
    p = 0
    while norm(g(λ + (δ / 2^p)))^2 > norm(g(λ))^2
      p += 1
    end

    λ += δ / 2^p
    inc = norm(δ)
    iter += 1
  end

  return (λ, Ẽ(λ), iter)
end

#=
println("---------- Testing ----------")

using Format

x = [0, 1, 2, 3, 4]
y = [3, 1, 0.5, 0.2, 0.05]

# Todo: check to see if this works, if x is a Vector{<:Real}
function f(x, coefs)
  a, b = coefs
  return a * exp(b * x)
end

h = 1e-7
λ = [1.0, 1.0]

λ, Ẽ, n = fit(f, λ, x, y, h)
a, b = λ[1], λ[2]

println("Gauss-Newton")
printfmtln("f(x) = {:.2f}⋅ℯ^({:.2f}⋅x)", a, b)
printfmtln("Ẽ = {:.2e}", Ẽ)
println("in $n iterations")
println()

λ = [2.0, 2.0]

λ, Ẽ, n = fit_damped(f, λ, x, y, h)
a, b = λ[1], λ[2]

println("Gauss-Newton damped")
printfmtln("f(x) = {:.2f}⋅ℯ^({:.2f}⋅x)", a, b)
printfmtln("Ẽ = {:.2e}", Ẽ)
println("in $n iterations")
=#
