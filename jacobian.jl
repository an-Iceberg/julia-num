include("partial.jl")

function D(f::Function, x⃗::Vector, h::Float64=1e-3)
  return stack([∂1_2(f, x⃗, i, h) for i in 1:length(x⃗)])
end

function D2(f::Function, x⃗::Vector, h::Float64=1e-3)
  return D(x -> D(f, x, h), x⃗, h)
end

function D3(f::Function, x⃗::Vector, h::Float64=1e-3)
  return D(x -> D2(f, x, h), x⃗, h)
end

function D4(f::Function, x⃗::Vector, h::Float64=1e-3)
  return D(x -> D3(f, x, h), x⃗, h)
end

# D¹

function D1_2(f::Function, x⃗::Vector{<:Real}, h::Float64=1e-3)::Matrix
  return stack([∂1_2(f, x⃗, i, h) for i in 1:length(x⃗)])
end

function D1_4(f::Function, x⃗::Vector{<:Real}, h::Float64=1e-3)::Matrix
  return stack([∂1_4(f, x⃗, i, h) for i in 1:length(x⃗)])
end

function D1_6(f::Function, x⃗::Vector{<:Real}, h::Float64=1e-3)::Matrix
  return stack([∂1_6(f, x⃗, i, h) for i in 1:length(x⃗)])
end

function D1_8(f::Function, x⃗::Vector{<:Real}, h::Float64=1e-3)::Matrix
  return stack([∂1_8(f, x⃗, i, h) for i in 1:length(x⃗)])
end

# D²

function D2_2(f::Function, x⃗::Vector{<:Real}, h::Float64=1e-3)::Matrix
  return stack([∂2_2(f, x⃗, i, h) for i in 1:length(x⃗)])
end

function D2_4(f::Function, x⃗::Vector{<:Real}, h::Float64=1e-3)::Matrix
  return stack([∂2_4(f, x⃗, i, h) for i in 1:length(x⃗)])
end

function D2_6(f::Function, x⃗::Vector{<:Real}, h::Float64=1e-3)::Matrix
  return stack([∂2_6(f, x⃗, i, h) for i in 1:length(x⃗)])
end

function D2_8(f::Function, x⃗::Vector{<:Real}, h::Float64=1e-3)::Matrix
  return stack([∂2_8(f, x⃗, i, h) for i in 1:length(x⃗)])
end

# D³

function D3_2(f::Function, x⃗::Vector{<:Real}, h::Float64=1e-3)::Matrix
  return stack([∂3_2(f, x⃗, i, h) for i in 1:length(x⃗)])
end

function D3_4(f::Function, x⃗::Vector{<:Real}, h::Float64=1e-3)::Matrix
  return stack([∂3_4(f, x⃗, i, h) for i in 1:length(x⃗)])
end

function D3_6(f::Function, x⃗::Vector{<:Real}, h::Float64=1e-3)::Matrix
  return stack([∂3_6(f, x⃗, i, h) for i in 1:length(x⃗)])
end

# D⁴

function D4_2(f::Function, x⃗::Vector{<:Real}, h::Float64=1e-3)::Matrix
  return stack([∂4_2(f, x⃗, i, h) for i in 1:length(x⃗)])
end

function D4_4(f::Function, x⃗::Vector{<:Real}, h::Float64=1e-3)::Matrix
  return stack([∂4_4(f, x⃗, i, h) for i in 1:length(x⃗)])
end

function D4_6(f::Function, x⃗::Vector{<:Real}, h::Float64=1e-3)::Matrix
  return stack([∂4_6(f, x⃗, i, h) for i in 1:length(x⃗)])
end

# D⁵

function D5_2(f::Function, x⃗::Vector{<:Real}, h::Float64=1e-3)::Matrix
  return stack([∂5_2(f, x⃗, i, h) for i in 1:length(x⃗)])
end

function D5_4(f::Function, x⃗::Vector{<:Real}, h::Float64=1e-3)::Matrix
  return stack([∂5_4(f, x⃗, i, h) for i in 1:length(x⃗)])
end

function D5_6(f::Function, x⃗::Vector{<:Real}, h::Float64=1e-3)::Matrix
  return stack([∂5_6(f, x⃗, i, h) for i in 1:length(x⃗)])
end

# D⁶

function D6_2(f::Function, x⃗::Vector{<:Real}, h::Float64=1e-3)::Matrix
  return stack([∂6_2(f, x⃗, i, h) for i in 1:length(x⃗)])
end

function D6_4(f::Function, x⃗::Vector{<:Real}, h::Float64=1e-3)::Matrix
  return stack([∂6_4(f, x⃗, i, h) for i in 1:length(x⃗)])
end

function D6_6(f::Function, x⃗::Vector{<:Real}, h::Float64=1e-3)::Matrix
  return stack([∂6_6(f, x⃗, i, h) for i in 1:length(x⃗)])
end

#=

# x = [1, 2, 3]
# y = copy(x)
# x[1] = 5
# display(x)
# display(y)

# function f⃗(x⃗)
#   x, y, z = x⃗[1], x⃗[2], x⃗[3]
#   return [
#     2x,
#     y,
#     0.5z^2
#   ]
# end

# function f(x)
#   x, y, z = x[1], x[2], x[3]
#   return [
#     2x,
#     y,
#     0.5z^2
#   ]
# end

# x⃗ = [1, 2, 3]
# display(f(x⃗))
# display(D1_2(f, x⃗))

# ---------

# function f1(input)
#   x, y, z = input
#   return x + y^2 - z^2 - 13
# end

# function f2(input)
#   x, y, z = input
#   return log(y / 4) + ℯ^(0.5z - 1) - 1
# end

# function f3(input)
#   x, y, z = input
#   return (y - 3)^2 - z^3 + 7
# end

# function f(λ)
#   return [f1(λ), f2(λ), f3(λ)]
# end

# f1(x) = x[1] + x[2]^2 - x[3]^2 - 13
# f2(x) = log(x[2] / 4) + ℯ^(0.5x[3] - 1) - 1
# f3(x) = (x[2] - 3)^2 - x[3]^3 + 7

# f(λ) = [f1(λ), f2(λ), f3(λ)]
# x⃗ = [1.5, 3, 2.5]

# display(D1_2(f, x⃗))
# display(D1_4(f, x⃗))
# display(D1_6(f, x⃗))
# display(D1_8(f, x⃗))

# Todo: use Format

using Printf

@printf "--------------------- Jacobian (Df) ---------------------\n\n"

f1(x) = x[1] + x[2]^2 - x[3]^2 - 13
f2(x) = log(x[2] / 4) + ℯ^(0.5x[3] - 1) - 1
f3(x) = (x[2] - 3)^2 - x[3]^3 + 7

f(λ) = [f1(λ), f2(λ), f3(λ)]
x⃗ = [1.5, 3, 2.5]

display(D1_2(f, x⃗))
display(D1_4(f, x⃗))
display(D1_6(f, x⃗))
display(D1_8(f, x⃗))
=#
