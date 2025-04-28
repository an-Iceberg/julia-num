using Format

include("derivative.jl")
include("integral.jl")

f(x) = cos(x)
a, b, h = 1, 5, 0.001

printfmtln("Fläche: {:.4f}", ∫_2(a, b, f, h))
printfmtln("Rotationsvolumen: {:.4f}", π * ∫_2(a, b, x -> f(x)^2, h))
printfmtln("Bogenlänge: {:.4f}", ∫_2(a, b, x -> sqrt(1 + d1_4(f, x, h)^2), h))
