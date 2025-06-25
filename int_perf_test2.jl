include("src/integral.jl")

using Dates

f(x) = log(x)
a = 0.1
b = 5_000_000

F(x) = x * (log(abs(x)) - 1)
I() = F(b) - F(a)

t1 = now()
result = ∫_2(a, b, f)
t2 = now()
println("∫_2: $(abs(t1 - t2) / Millisecond(1_000))  ε = $(round(abs(I() - result), sigdigits=2))")

t1 = now()
result = ∫_3(a, b, f)
t2 = now()
println("∫_3: $(abs(t1 - t2) / Millisecond(1_000))  ε = $(round(abs(I() - result), sigdigits=2))")

t1 = now()
result = ∫_4(a, b, f)
t2 = now()
println("∫_4: $(abs(t1 - t2) / Millisecond(1_000))  ε = $(round(abs(I() - result), sigdigits=2))")

t1 = now()
result = ∫_6(a, b, f)
t2 = now()
println("∫_6: $(abs(t1 - t2) / Millisecond(1_000))  ε = $(round(abs(I() - result), sigdigits=2))")
