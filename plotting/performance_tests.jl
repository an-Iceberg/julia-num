include("functions.jl")

using DataFrames, Dates, CSV, Statistics

# Todo: do statistics

iter_counts = [1_000, 2_000, 5_000, 10_000, 20_000, 50_000, 100_000, 200_000, 500_000, 1_000_000, 2_000_000, 5_000_000]
slow_times = []
fast_times = []

for _ in 1:10
  ∂1_4_slow(f, x⃗, 1)
  ∂1_4_fast(f, x⃗, 1)
end

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

println("Saving data to file…")
data = DataFrame(test_size=iter_counts, slow_times=slow_times, fast_times=fast_times)
CSV.write("plotting/data.csv", data)
println("Done…")
