# Integrals

"""
Calculates the integral of `f` between `a` and `b` using Simpson's 1/3 rule, so
using a polynomial of 2ⁿᵈ degree.
"""
function ∫_2(a, b, f, h)
  n = ceil((b - a) / h)
  x(i) = a + (i * h)

  Σ₁ = sum([f(x(2i - 1)) for i in 1:(n/2)])
  Σ₂ = sum([f(x(2i)) for i in 1:((n/2)-1)])

  return (1 / 3) * h * (f(a) + 4Σ₁ + 2Σ₂ + f(b))
end

"""
Calculates the integral of `f` between `a` and `b` using Simpson's 3/8 rule, so
using a polynomial of 3ʳᵈ degree.
"""
function ∫_3(a, b, f, h)
  n = ceil((b - a) / h)
  x(i) = a + (i * h)

  Σ₁ = collect(1:n)
  filter!(i -> i % 3 != 0, Σ₁)
  Σ₁ = sum([f(x(i)) for i in Σ₁])
  Σ₂ = sum([f(x(3i)) for i in 1:((n/3)-1)])

  return (3 / 8) * h * (f(a) + 3Σ₁ + 2Σ₂ + f(b))

  # Σ = f(a) + f(b)
  # for i in 1:n
  #   if i % 3 == 0
  #     Σ += 2 * f(x(i))
  #   else
  #     Σ += 3 * f(x(i))
  #   end
  # end
  # return (3 / 8) * h * Σ

  # return (3 / 8) * h *
  #        sum([f(x(3i - 3)) + f(x(3i - 2)) + f(x(3i - 1)) + f(x(3i)) for i in 1:(n / 3 - 1)])
end

"""
Calculates the integral of `f` between `a` and `b` using Boole's rule, so
using a polynomial of 4ʳᵈ degree.
"""
function ∫_4(a, b, f, h)
  n = ceil((b - a) / h)
  x(i) = a + (i * h)

  Σ₁ = sum([f(x(1 + 2i)) for i in 0:((n/2)-1)])
  Σ₂ = sum([f(x(2 + 4i)) for i in 0:((n/4)-1)])
  Σ₃ = sum([f(x(4 + 4i)) for i in 0:((n/4)-2)])

  return (2 / 45) * h * (7 * (f(a) + f(b)) + 32Σ₁ + 12Σ₂ + 14Σ₃)
end
