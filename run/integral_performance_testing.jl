# Todo: proper performance analysis, meaning 100 samples for each iteration size and and box plots

function ∫_3_list(a::Real, b::Real, f::Function, h::Real=1e-2)::Real
  n = ceil((b - a) / h)
  x(i::Real)::Real = a + (i * h)
  return (3 / 8) * h * sum([f(x(3i - 3)) + 3f(x(3i - 2)) + 3f(x(3i - 1)) + f(x(3i)) for i in 1:(n/3)])
end

function ∫_3_loop(a::Real, b::Real, f::Function, h::Real=1e-2)::Real
  n = ceil((b - a) / h)
  x(i::Real)::Real = a + (i * h)

  Σ = f(a) + f(b)
  for i in 1:n
    if mod(i, 3) == 0
      Σ += 2 * f(x(i))
    else
      Σ += 3 * f(x(i))
    end
  end

  return (3 / 8) * h * Σ
end

# Terrible performance 😆
function ∫_3_map(a::Real, b::Real, f::Function, h::Real=1e-2)::Real
  n = ceil((b - a) / h)
  x(i::Real)::Real = a + (i * h)

  Σ₁ = collect(1:n)
  filter!(i -> i % 3 != 0, Σ₁)
  # filter!(i -> mod(i, 3) != 0, Σ₁)
  Σ₁ = sum([f(x(i)) for i in Σ₁])
  Σ₂ = sum([f(x(3i)) for i in 1:((n/3)-1)])

  return (3 / 8) * h * (f(a) + 3Σ₁ + 2Σ₂ + f(b))
end

f(x) = cos(x)
a = 0.0
b = 5.0

@time "single ∫_list" ∫_3_list(a, b, f)
@time "single ∫_loop" ∫_3_loop(a, b, f)
@time "single ∫_map" ∫_3_map(a, b, f)
println()

iter_count = 750_000

@time "loop ∫_list" begin
  for _ in 1:iter_count
    ∫_3_list(a, b, f)
  end
end

@time "loop ∫_loop" begin
  for _ in 1:iter_count
    ∫_3_loop(a, b, f)
  end
end

@time "loop ∫_map" begin
  for _ in 1:iter_count
    ∫_3_map(a, b, f)
  end
end

println()

using Dates

t1 = now()
for _ in 1:iter_count
  ∫_3_list(a, b, f)
end
t2 = now()
Δt = abs(t1 - t2)
println("∫_3_list: $Δt")

t1 = now()
for _ in 1:iter_count
  ∫_3_loop(a, b, f)
end
t2 = now()
Δt = abs(t1 - t2)
println("∫_3_loop: $Δt")

t1 = now()
for _ in 1:iter_count
  ∫_3_map(a, b, f)
end
t2 = now()
Δt = abs(t1 - t2)
println("∫_3_map: $Δt")
