
# Todo: proper performance analysis, meaning 100 samples for each iteration size and and box plots

function ∂1_2_slow(f::Function, x⃗::Vector{<:Real}, i::Int, h::Real=1e-3)::Union{Real,Vector{<:Real}}
  h⃗ = zeros(length(x⃗))
  h⃗[i] = h
  return (-0.5f(x⃗ - h⃗) + 0.5f(x⃗ + h⃗)) / h
end

function ∂1_2_fast(f::Function, x⃗::Vector{<:Real}, i::Int, h::Real=1e-3)::Union{Real,Vector{<:Real}}
  # This is for performance reasons. Each time we need to evalute f we only modify the relevant
  # vector entry instead of creating another vector.
  x_local = copy(x⃗)

  # Info: this is a closure and captures its environment. That's why it can't be refactored to be
  # placed outside of this function's scope
  function x(h)
    x_local[i] = x⃗[i]
    x_local[i] += h
    return x_local
  end

  return (-0.5f(x(-h)) + 0.5f(x(h))) / h
end

function ∂1_4_slow(f::Function, x⃗::Vector{<:Real}, i::Int, h::Real=1e-3)::Union{Real,Vector{<:Real}}
  h⃗ = zeros(length(x⃗))
  h⃗[i] = h
  return ((1 / 12)f(x⃗ - 2h⃗) - (2 / 3)f(x⃗ - h⃗) + (2 / 3)f(x⃗ + h⃗) - (1 / 12)f(x⃗ + 2h⃗)) / h
end

function ∂1_4_fast(f::Function, x⃗::Vector{<:Real}, i::Int, h::Real=1e-3)::Union{Real,Vector{<:Real}}
  x_local = copy(x⃗)

  # Info: this is a closure and captures its environment. That's why it can't be refactored to be
  # placed outside of this function's scope
  function x(h)
    x_local[i] = x⃗[i]
    x_local[i] += h
    return x_local
  end

  return ((1 / 12)f(x(-2h)) - (2 / 3)f(x(-h)) + (2 / 3)f(x(h)) - (1 / 12)f(x(2h))) / h
end

x⃗ = [Float64(π), Float64(π)]
f(x⃗) = sin(x⃗[1]) + cos(x⃗[2])

# @time "single ∂1_2_slow" begin
#   ∂1_2_slow(f, x⃗, 1)
# end
# @time "single ∂1_2_fast" begin
#   ∂1_2_fast(f, x⃗, 1)
# end

# @time "loop ∂1_2_slow" begin
#   for _ in 1:1000
#     ∂1_2_slow(f, x⃗, 1)
#   end
# end

# @time "loop ∂1_2_fast" begin
#   for _ in 1:1000
#     ∂1_2_fast(f, x⃗, 1)
#   end
# end

# println()

# @time "single ∂1_4_slow" begin
#   ∂1_4_slow(f, x⃗, 1)
# end
# @time "single ∂1_4_fast" begin
#   ∂1_4_fast(f, x⃗, 1)
# end

# @time "loop ∂1_4_slow" begin
#   for _ in 1:1000
#     ∂1_4_slow(f, x⃗, 1)
#   end
# end

# @time "loop ∂1_4_fast" begin
#   for _ in 1:1000
#     ∂1_4_fast(f, x⃗, 1)
#   end
# end

# @showtime
# @timev
# @timed
# @elapsed

println()

using Dates

# iter_count = 5_000_000

# t1 = now()
# for _ in 1:iter_count
#   ∂1_2_slow(f, x⃗, 1)
# end
# t2 = now()
# Δt = abs(t1 - t2) / Millisecond(1000)
# println("∂1_2_slow: $Δt")

# t1 = now()
# for _ in 1:iter_count
#   ∂1_2_fast(f, x⃗, 1)
# end
# t2 = now()
# Δt = abs(t1 - t2) / Millisecond(1000)
# println("∂1_2_fast: $Δt")

# t1 = now()
# for _ in 1:iter_count
#   ∂1_4_slow(f, x⃗, 1)
# end
# t2 = now()
# Δt = abs(t1 - t2) / Millisecond(1000)
# println("∂1_4_slow: $Δt")

# t1 = now()
# for _ in 1:iter_count
#   ∂1_4_fast(f, x⃗, 1)
# end
# t2 = now()
# Δt = abs(t1 - t2) / Millisecond(1000)
# println("∂1_4_fast: $Δt")

# %%%%%%%%%%%%%%%%%%%% Stats %%%%%%%%%%%%%%%%%%%%

using CairoMakie
using DataFrames
using Format

# https://docs.makie.org/dev/reference/plots/barplot
# https://discourse.julialang.org/t/dodged-barplot-got-thin-out-in-makie/87923/5
# https://discourse.julialang.org/t/need-help-in-formatting-bar-plot-in-makie/80496

# Todo: use Dataframes
# Todo: test statistically for each iter count using 100 samples
iter_counts = [100, 200, 500, 1_000, 2_000, 5_000, 10_000, 20_000, 50_000, 100_000, 200_000, 500_000]
slow_times = []
fast_times = []

println("Running performances tests…")
for iter_count in iter_counts
  t1 = now()
  for _ in 1:iter_count
    ∂1_4_slow(f, x⃗, 1)
  end
  t2 = now()
  Δt_slow = abs(t1 - t2) / Millisecond(1000)
  append!(slow_times, Δt_slow)

  t1 = now()
  for _ in 1:iter_count
    ∂1_4_fast(f, x⃗, 1)
  end
  t2 = now()
  Δt_fast = abs(t1 - t2) / Millisecond(1000)
  append!(fast_times, Δt_fast)
end
println("Done…")
println()

println("Creating plot…")

x_labels = []

for n in iter_counts
  str = ""

  if n >= 1_000 && n < 1_000_000
    str = format("{:d}k", n / 1_000)
  elseif n >= 1_000_000 && n < 1_000_000_000
    str = format("{:d}M", n / 1_000_000)
  else
    str = format("{}", n)
  end

  push!(x_labels, str)
end

figure = Figure()
ax = Axis(
  figure[1, 1],
  title="Benchmark",
  xlabel="# of function calls",
  ylabel="duration in seconds",
  xticks=(1:length(iter_counts), x_labels), #=[replace(cfmt("%'d", n), "," => "'") for n in iter_counts]=#
  yscale=log
)
categories = repeat(1:length(iter_counts), 2)
groups = vcat(repeat([1], length(iter_counts)), repeat([2], length(iter_counts)))
# println(categories)
# println(groups)
barplot!(
  ax,
  categories,
  vcat(slow_times, fast_times),
  dodge=groups,
  color=Makie.wong_colors()[groups],
  # strokewidth=1,
  # width=0.5,
  # gap=0
)
# axislegend(position=:lt)
# Legend(figure[1, 2], [PolyElement(polycolor=Makie.wong_colors()[i]) for i in 1:2], ["slow impl.", "fast impl."], "Groups")
save("benchmark.png", figure)
println("Done…")

# Todo: research this
# using Format
