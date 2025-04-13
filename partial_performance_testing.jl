
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

@time "single ∂1_2_slow" begin
  ∂1_2_slow(f, x⃗, 1)
end
@time "single ∂1_2_fast" begin
  ∂1_2_fast(f, x⃗, 1)
end

@time "loop ∂1_2_slow" begin
  for _ in 1:1000
    ∂1_2_slow(f, x⃗, 1)
  end
end

@time "loop ∂1_2_fast" begin
  for _ in 1:1000
    ∂1_2_fast(f, x⃗, 1)
  end
end

println()

@time "single ∂1_4_slow" begin
  ∂1_4_slow(f, x⃗, 1)
end
@time "single ∂1_4_fast" begin
  ∂1_4_fast(f, x⃗, 1)
end

@time "loop ∂1_4_slow" begin
  for _ in 1:1000
    ∂1_4_slow(f, x⃗, 1)
  end
end

@time "loop ∂1_4_fast" begin
  for _ in 1:1000
    ∂1_4_fast(f, x⃗, 1)
  end
end

# @showtime
# @timev
# @timed
# @elapsed

println()

using Dates

iter_count = 5_000_000

t1 = now()
for _ in 1:iter_count
  ∂1_2_slow(f, x⃗, 1)
end
t2 = now()
Δt = abs(t1 - t2)
println("∂1_2_slow: $Δt")

t1 = now()
for _ in 1:iter_count
  ∂1_2_fast(f, x⃗, 1)
end
t2 = now()
Δt = abs(t1 - t2)
println("∂1_2_fast: $Δt")

t1 = now()
for _ in 1:iter_count
  ∂1_4_slow(f, x⃗, 1)
end
t2 = now()
Δt = abs(t1 - t2)
println("∂1_4_slow: $Δt")

t1 = now()
for _ in 1:iter_count
  ∂1_4_fast(f, x⃗, 1)
end
t2 = now()
Δt = abs(t1 - t2)
println("∂1_4_fast: $Δt")
