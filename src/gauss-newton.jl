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
function fit(f::Function, λ0, x, y, h::Real=1e-3, n::Int=100)
  D = ForwardDiff.jacobian
  g(λ) = [y - f(x, λ) for (x, y) in zip(x, y)]
  Ẽ(λ) = norm(g(λ))^2

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
function fit_damped(f::Function, λ₀, x, y, h::Real=1e-3, n::Integer=100)
  D = ForwardDiff.jacobian
  g(λ) = [y - f(x, λ) for (x, y) in zip(x, y)]
  Ẽ(λ) = norm(g(λ))^2

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
using Printf

x = [0, 1, 2, 3, 4]
y = [3, 1, 0.5, 0.2, 0.05]

function f(x, coefs)
  a, b = coefs
  return a * exp(b * x)
end

h = 1e-7
λ = [1.0, 1.0]

λ, Ẽ, n = fit(f, λ, x, y, h)
a, b = λ[1], λ[2]

println("Gauss-Newton")
@printf("f(x) = %.2f⋅ℯ^%.2f⋅x\n", a, b)
@printf("Ẽ = %.2f\n", Ẽ)
println("in $n iterations")
println()

λ = [2.0, 2.0]

λ, Ẽ, n = fit_damped(f, λ, x, y, h)
a, b = λ[1], λ[2]

println("Gauss-Newton damped")
@printf("f(x) = %.2f⋅ℯ^%.2f⋅x\n", a, b)
@printf("Ẽ = %.2f\n", Ẽ)
println("in $n iterations")
=#

#=
using Printf

x = [[0, 0], [3, 1], [2, 2], [1, 3], [3, 4]]
y = [3, 1, 0.5, 0.2, 0.05]

function fn(x, coefs)
  x, y = x
  a, b, c = coefs
  return a * x + b * y + c
end

a, b, c = 0, 0, 0
params = [a, b, c]

params, error, iter_count = fit_damped(fn, params, x, y, 1e-7)
a, b, c = params

prec = 4
@printf("f(x, y) = %.*fx + %.*fy + %.*f\n", prec, a, prec, b, prec, c)
@printf("Ẽ = %.2e\n", error)
println("in $iter_count iterations")
=#

#=
# Todo: function that returns a vector
using Printf

x = [[0, 0], [3, 1], [2, 2], [1, 3], [3, 4]]
y = [3, 1, 0.5, 0.2, 0.05]

function fn(x, coefs)
  x, y = x
  a, b, c = coefs
  return a * x + b * y + c
end

a, b, c = 0, 0, 0
params = [a, b, c]

params, error, iter_count = fit_damped(fn, params, x, y, 1e-7)
a, b, c = params

prec = 4
@printf("f(x, y) = %.*fx + %.*fy + %.*f\n", prec, a, prec, b, prec, c)
@printf("Ẽ = %.2e\n", error)
println("in $iter_count iterations")
=#
