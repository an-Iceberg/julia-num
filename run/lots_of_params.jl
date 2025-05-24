using Printf

include("../src/gauss-newton_fin_diff.jl")

x = [1, 2, 3, 4, 5]
y = [1, 1.5, 0.5, 0.75, 1]

function f(x, coefs)
  a, b, c, d, e, f, g, h, i, j, k, l = coefs

  sin_part(a, b, c, d) = a * sin(b * x + c) + d

  return sin_part(a, b, c, d) + sin_part(e, f, g, h) + sin_part(i, j, k, l)
end

λ = fill(1, 12)

λ, Ẽ, n = fit_damped(f, λ, x, y)

using CairoMakie

x_line = collect(LinRange(0, 6, 100))
y_line = [f(x, λ) for x in x_line]

figure = Figure()
ax = Axis(figure[1, 1]; title="Data fit")
lines!(ax, x_line, y_line; label="f(x, λ)")
scatter!(ax, x, y; label="data", color=:orange)
axislegend(; position=:rb)

x_line = collect(LinRange(-60, 65, 600))
y_line = [f(x, λ) for x in x_line]

ax2 = Axis(figure[2, 1]; title="The bigger picture")
lines!(ax2, x_line, y_line)
scatter!(ax2, x, y; color=:orange)

save("lots_of_params.png", figure)
