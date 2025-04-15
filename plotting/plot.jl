using DataFrames, CSV, CairoMakie, Format

data = DataFrame(CSV.File("plotting/data.csv"))

test_size = data[!, "test_size"]
slow_times = data[!, "slow_times"]
fast_times = data[!, "fast_times"]
x_labels = []

for n in test_size
  str = ""

  if n >= 1_000 && n < 1_000_000
    str = format("{:d}k", n / 1_000)
  elseif n >= 1_000_000 && n < 1_000_000_000
    str = format("{:d}M", n / 1_000_000)
  elseif n >= 1_000_000_000 && n < 1_000_000_000_000
    str = format("{:d}G", n / 1_000_000_000)
  else
    str = format("{}", n)
  end

  push!(x_labels, str)
end

figure = Figure(title="Performance benchmark")

ax = Axis(
  figure[1, 1],
  title="Linear scale",
  xlabel="# of function calls",
  ylabel="duration (s)",
  xticks=(1:length(test_size), x_labels),
)
categories = repeat(1:length(test_size), 2)
groups = vcat(repeat([1], length(test_size)), repeat([2], length(test_size)))
barplot!(
  ax,
  categories,
  vcat(slow_times, fast_times),
  dodge=groups,
  color=Makie.wong_colors()[groups],
)

ax_log = Axis(
  figure[2, 1],
  title="Logarithmic scale",
  xlabel="# of function calls",
  ylabel="duration (s)",
  xticks=(1:length(test_size), x_labels),
  yscale=log
)
categories = repeat(1:length(test_size), 2)
groups = vcat(repeat([1], length(test_size)), repeat([2], length(test_size)))
barplot!(
  ax_log,
  categories,
  vcat(slow_times, fast_times),
  dodge=groups,
  color=Makie.wong_colors()[groups],
)

save("plotting/benchmark.png", figure)
