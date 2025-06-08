# Todo: add links to wikipedia formulas

"""
Calculates the integral of `f` between `a` and `b` using
[Simpson's ⅓ rule](https://en.wikipedia.org/wiki/Simpson%27s_rule#Composite_Simpson's_1/3_rule)
, so using polynomials of 2ⁿᵈ degree.
"""
function ∫_2(a::Real, b::Real, f::Function, h::Real=1e-2)::Real
  # Todo: don't ceil and don't ::Int. Do that at the array comprehension
  n = ceil((b - a) / h)
  x(i) = a + (i * h)

  # Todo: performance test, which one of these 2 is faster. They each give same precision

  # Σ₁ = sum([f(x(2i - 1)) for i in 1:(n/2)])
  # Σ₂ = sum([f(x(2i)) for i in 1:((n/2)-1)])

  # return (1 / 3) * h * (f(a) + 4Σ₁ + 2Σ₂ + f(b))

  #         ⅓      ‧       (1      +       4         +     1       )     =        2
  return (1 / 3) * h * sum([f(x(2i - 2)) + 4f(x(2i - 1)) + f(x(2i)) for i in 1:(n / 2)])
end

"""
Calculates the integral of `f` between `a` and `b` using
[Simpson's ⅜ rule](https://en.wikipedia.org/wiki/Simpson%27s_rule#Composite_Simpson's_3/8_rule)
, so using polynomials of 3ʳᵈ degree.
"""
function ∫_3(a::Real, b::Real, f::Function, h::Real=1e-2)::Real
  n = ceil((b - a) / h)
  x(i) = a + (i * h)

  # modulus = n % 3
  # if modulus == 1
  #   n -= 1
  # elseif modulus == 2
  #   n += 1
  # end

  # Note: none of these work 😭

  # Σ₁ = collect(1:n)
  # filter!(i -> i % 3 != 0, Σ₁)
  # # filter!(i -> mod(i, 3) != 0, Σ₁)
  # Σ₁ = sum([f(x(i)) for i in Σ₁])
  # Σ₂ = sum([f(x(3i)) for i in 1:((n/3)-1)])

  # return (3 / 8) * h * (f(a) + 3Σ₁ + 2Σ₂ + f(b))

  # Σ = f(a) + f(b)
  # for i in 1:n
  #   if mod(i, 3) == 0
  #     Σ += 2 * f(x(i))
  #   else
  #     Σ += 3 * f(x(i))
  #   end
  # end
  # return (3 / 8) * h * Σ

  #! format: off
  return (3/8) * h * sum([f(x(3i-3)) + 3f(x(3i-2)) + 3f(x(3i-1)) + f(x(3i)) for i in 1:(n/3)])
  #! format: on

  # return (3 / 8) * h * sum([f(x(3i - 3)) + f(x(3i - 2)) + f(x(3i - 1)) + f(x(3i)) for i in 1:(n/3-1)])
end

"""
Calculates the integral of `f` between `a` and `b` using
[Boole's rule](https://en.wikipedia.org/wiki/Finite_difference_coefficient)
, so using polynomials of 4ᵗʰ degree.
"""
function ∫_4(a::Real, b::Real, f::Function, h::Real=1e-2)::Real
  n = ceil((b - a) / h)
  x(i) = a + (i * h)

  # Todo: performance test, which one of these 2 is faster. They each give same precision

  # Σ₁ = sum([f(x(1 + 2i)) for i in 0:((n/2)-1)])
  # Σ₂ = sum([f(x(2 + 4i)) for i in 0:((n/4)-1)])
  # Σ₃ = sum([f(x(4 + 4i)) for i in 0:((n/4)-2)])

  # return (2 / 45) * h * (7 * (f(a) + f(b)) + 32Σ₁ + 12Σ₂ + 14Σ₃)

  #! format: off
  return (2/45) * h * sum([7f(x(4i-4)) + 32f(x(4i-3)) + 12f(x(4i-2)) + 32f(x(4i-1)) + 7f(x(4i)) for i in 1:(n/4)])
  #! format: on
end

"""
Calculates the integral of `f` between `a` and `b` using
[Weddle's rule]()
, so using polynomials of 6ᵗʰ degree.
"""
function ∫_6(a::Real, b::Real, f::Function, h::Real=1e-2)::Real
  n = ceil((b - a) / h)
  x(i) = a + (i * h)

  #! format: off
  return (3/10) * h * sum([f(x(6i-6)) + 5f(x(6i-5)) + f(x(6i-4)) + 6f(x(6i-3)) + f(x(6i-2)) + 5f(x(6i-1)) + f(x(6i)) for i in 1:(n/6)])
  #! format: on
end

#=
using Printf

f(x) = cos(x)
F(x) = sin(x)
a = 0.0
b = 5.0

@printf("Simpson's 1/3: ε = %.2e\n", abs(∫_2(a, b, f) - (F(b) - F(a))))
@printf("Simpson's 3/8: ε = %.2e\n", abs(∫_3(a, b, f) - (F(b) - F(a))))
@printf("Boole's      : ε = %.2e\n", abs(∫_4(a, b, f) - (F(b) - F(a))))
@printf("Weddle's     : ε = %.2e\n", abs(∫_6(a, b, f) - (F(b) - F(a))))
=#
