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

(Fitted coefficients, mean square error, \\# of iterations)
"""
function fit(
  f::Function,
  λ₀::Vector{<:Real},
  x::Vector{<:Real},
  y::Vector{<:Real},
  h::Real=1e-3,
  n::Int=200,
)::Tuple{Vector{<:Real}, Real, Int}
  g(λ::Vector{<:Real})::Vector{<:Real} = [y - f(x, λ) for (x, y) in zip(x, y)]
  Ẽ(λ::Vector{<:Real})::Real = norm(g(λ))^2

  increment = h + 1
  iter_count = 0
  λ = copy(λ₀)

  while iter_count < n && increment > h
    Q, R = qr(D1_4(g, λ, h))
    δ = R \ (.-transpose(Matrix(Q)) * g(λ))
    λ += δ
    increment = norm(δ)
    iter_count += 1
  end

  return (λ, Ẽ(λ), iter_count)
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
  λ₀::Vector{<:Real},
  x::Vector{<:Real},
  y::Vector{<:Real},
  h::Real=1e-3,
  n::Int=200,
)::Tuple{Vector{<:Real}, Real, Int}
  g(λ::Vector{<:Real})::Vector{<:Real} = [y - f(x, λ) for (x, y) in zip(x, y)]
  Ẽ(λ::Vector{<:Real})::Real = norm(g(λ))^2

  increment = h + 1
  iter_count = 0
  λ = copy(λ₀)

  while iter_count < n && increment > h
    Q, R = qr(D1_4(g, λ, h))
    δ = R \ (.-transpose(Matrix(Q)) * g(λ))

    p = 0
    while norm(g(λ + (δ / 2^p)))^2 > norm(g(λ))^2
      p += 1
    end

    λ += δ / 2^p
    increment = norm(δ)
    iter_count += 1
  end

  return (λ, Ẽ(λ), iter_count)
end

#=
using Format

x = [0, 1, 2, 3, 4]
y = [3, 1, 0.5, 0.2, 0.05]

function f(x::Real, coefs::Vector{<:Real})::Real
  a, b = coefs
  return a * ℯ^(b * x)
end

h = 1e-7
λ = [1.0, 1.0]

λ, Ẽ, n = fit(f, λ, x, y, h)
a, b = λ[1], λ[2]

println("Gauss-Newton")
printfmtln("f(x) = {:.2f}⋅ℯ^({:.2f}⋅x)", a, b)
printfmtln("Ẽ = {:.2e}", Ẽ)
println("in $n iterations")

x_line = collect(LinRange(x[1], x[end], 200))
y_line = [f(x, λ) for x in x_line]

using CairoMakie

figure = Figure()

ax = Axis(figure[1, 1]; title="Gauss-Newton undamped example")
lines!(ax, x_line, y_line; label="fitted function")
scatter!(ax, x, y; label="data", color=:orange)

axislegend(; position=:rt)

save("gauss_newton_undamped.png", figure)

println()

λ = [2.0, 2.0]

λ, Ẽ, n = fit_damped(f, λ, x, y, h)
a, b = λ[1], λ[2]

println("Gauss-Newton damped")
printfmtln("f(x) = {:.2f}⋅ℯ^({:.2f}⋅x)", a, b)
printfmtln("Ẽ = {:.2e}", Ẽ)
println("in $n iterations")

x_line = collect(LinRange(x[1], x[end], 200))
y_line = [f(x, λ) for x in x_line]

figure = Figure()

ax = Axis(figure[1, 1]; title="Gauss-Newton damped example")
lines!(ax, x_line, y_line; label="fitted function")
scatter!(ax, x, y; label="data", color=:orange)

axislegend(; position=:rt)

save("gauss_newton_damped.png", figure)
=#
