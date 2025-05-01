include("../derivative.jl")

using Format

println("--------------------- Differentiation (df) ---------------------")
println()

x = Float64(π)
f(x) = exp(x)
d1f(x) = exp(x)
d2f(x) = exp(x)
d3f(x) = exp(x)
d4f(x) = exp(x)
d5f(x) = exp(x)
d6f(x) = exp(x)

printfmtln("d1_2: {:.2e}", abs(d1f(x) - d1_2(f, x)))
printfmtln("d1_4: {:.2e}", abs(d1f(x) - d1_4(f, x)))
printfmtln("d1_6: {:.2e}", abs(d1f(x) - d1_6(f, x)))
printfmtln("d1_8: {:.2e}", abs(d1f(x) - d1_8(f, x)))
println()
printfmtln("d2_2: {:.2e}", abs(d2f(x) - d2_2(f, x)))
printfmtln("d2_4: {:.2e}", abs(d2f(x) - d2_4(f, x)))
printfmtln("d2_6: {:.2e}", abs(d2f(x) - d2_6(f, x)))
printfmtln("d2_8: {:.2e}", abs(d2f(x) - d2_8(f, x)))
println()
printfmtln("d3_2: {:.2e}", abs(d3f(x) - d3_2(f, x)))
printfmtln("d3_4: {:.2e}", abs(d3f(x) - d3_4(f, x)))
printfmtln("d3_6: {:.2e}", abs(d3f(x) - d3_6(f, x)))
println()
printfmtln("d4_2: {:.2e}", abs(d4f(x) - d4_2(f, x, 0.01)))
printfmtln("d4_4: {:.2e}", abs(d4f(x) - d4_4(f, x, 0.01)))
printfmtln("d4_6: {:.2e}", abs(d4f(x) - d4_6(f, x, 0.01)))
println()
printfmtln("d5_2: {:.2e}", abs(d5f(x) - d5_2(f, x, 0.01)))
printfmtln("d5_4: {:.2e}", abs(d5f(x) - d5_4(f, x, 0.01)))
printfmtln("d5_6: {:.2e}", abs(d5f(x) - d5_6(f, x, 0.01)))
println()
printfmtln("d6_2: {:.2e}", abs(d6f(x) - d6_2(f, x, 1.0)))
printfmtln("d6_4: {:.2e}", abs(d6f(x) - d6_4(f, x, 1.0)))
printfmtln("d6_6: {:.2e}", abs(d6f(x) - d6_6(f, x, 0.1)))
println()

function d5(f::Function, x::Real, h::Real=1e-3)::Real
  return (-f(x - 3h) + 4f(x - 2h) - 5f(x - h) + 5f(x + h) - 4f(x + 2h) + f(x + 3h)) / 2h^5
end

function d6(f::Function, x::Real, h::Real=1e-3)::Real
  return (13f(x - 5h) - 190f(x - 4h) + 1_305f(x - 3h) - 4_680f(x - 2h) + 9_690f(x - h) - 12_276f(x) + 9_60f(x + h) - 4_680f(x + 2h) + 1_305f(x + 3h) - 190f(x + 4h) + 13f(x + 5h)) / 240h^6
end

printfmtln("d5: {:.2e}", abs(d5f(x) - d5(f, x, 0.1)))
printfmtln("d6: {:.2e}", abs(d6f(x) - d6(f, x, 0.1)))
