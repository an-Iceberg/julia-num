using Dates

function ∫_2_1(a::Real, b::Real, f::Function, h::Real=1e-2)::Real
  n = ceil((b - a) / h)
  x(i) = a + (i * h)
  return (1 / 3) * h * sum([f(x(2i - 2)) + 4f(x(2i - 1)) + f(x(2i)) for i in 1:(n / 2)])
end

function ∫_2_2(a::Real, b::Real, f::Function, h::Real=1e-2)::Real
  n = ceil((b - a) / h)
  x(i) = a + (i * h)
  Σ₁ = sum([f(x(2i - 1)) for i in 1:(n / 2)])
  Σ₂ = sum([f(x(2i)) for i in 1:((n / 2) - 1)])
  return (1 / 3) * h * (f(a) + 4Σ₁ + 2Σ₂ + f(b))
end

t1 = now()
∫_2_1(0.1, 1_000_000, log)
t2 = now()
println("array comprehension: $(abs(t1 - t2) / Millisecond(1_000))")

t1 = now()
∫_2_2(0.1, 1_000_000, log)
t2 = now()
println("partial sums: $(abs(t1 - t2) / Millisecond(1_000))")
