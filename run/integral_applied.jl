using Format

include("../src/derivative.jl")
include("../src/integral.jl")

begin
  f(x) = cos(x)
  a, b, h = 1, 5, 0.001

  printfmtln("Fläche: {:.4f}", ∫_2(a, b, f, h))
  printfmtln("Rotationsvolumen: {:.4f}", π * ∫_2(a, b, x -> f(x)^2, h))
  printfmtln("Bogenlänge: {:.4f}", ∫_2(a, b, x -> sqrt(1 + d1_4(f, x, h)^2), h))
end

println()

begin
  f(x) = -x^2 + 4x + 3
  g(x) = -x^3 + 7x^2 - 10x + 5
  a, b = 1, 2

  printfmtln("{:.4f}", ∫_2(a, b, x -> abs(f(x) - g(x))))
  printfmtln("{:.4f}", 49 / 12)
end

println()

begin
  f(x) = sqrt(x)
  g(x) = sqrt(x + 1)
  a, b = 0, 4

  printfmtln("{:.4f}", ∫_2(a, b, x -> abs(f(x) - g(x))))
  printfmtln("{:.4f}", 10sqrt(5) / 3 - 6)
end

println()

begin
  f(x) = sin(x) * cos(x)
  g(x) = sin(x)
  a, b = 0, π

  printfmtln("{:.4f}", ∫_2(a, b, x -> abs(f(x) - g(x))))
  printfmtln("{:.4f}", 2)
end
