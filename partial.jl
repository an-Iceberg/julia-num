include("deriv_rework.jl")

# ∂¹

function ∂1_2(f::Function, x⃗::Vector{<:Real}, i::Int, h::Float64=1e-3)
  h⃗ = zeros(length(x⃗))
  h⃗[i] = h
  return (-(1 / 2)f(x⃗ - h⃗) + (1 / 2)f(x⃗ + h⃗)) / h
end

function ∂1_4(f::Function, x⃗::Vector{<:Real}, i::Int, h::Float64=1e-3)
  h⃗ = zeros(length(x⃗))
  h⃗[i] = h
  return ((1 / 12)f(x⃗ - 2h⃗) - (2 / 3)f(x⃗ - h⃗) + (2 / 3)f(x⃗ + h⃗) - (1 / 12)f(x⃗ + 2h⃗)) / h
end

function ∂1_6(f::Function, x⃗::Vector{<:Real}, i::Int, h::Float64=1e-3)
  h⃗ = zeros(length(x⃗))
  h⃗[i] = h
  return (-(1 / 60)f(x⃗ - 3h⃗) + (3 / 20)f(x⃗ - 2h⃗) - (3 / 4)f(x⃗ - h⃗) + (3 / 4)f(x⃗ + h⃗) - (3 / 20)f(x⃗ + 2h⃗) + (1 / 60)f(x⃗ + 3h⃗)) / h
end

function ∂1_8(f::Function, x⃗::Vector{<:Real}, i::Int, h::Float64=1e-3)
  h⃗ = zeros(length(x⃗))
  h⃗[i] = h
  return ((1 / 280)f(x⃗ - 4h⃗) - (4 / 105)f(x⃗ - 3h⃗) + (1 / 5)f(x⃗ - 2h⃗) - (4 / 5)f(x⃗ - h⃗) + (4 / 5)f(x⃗ + h⃗) - (1 / 5)f(x⃗ + 2h⃗) + (4 / 105)f(x⃗ + 3h⃗) - (1 / 280)f(x⃗ + 4h⃗)) / h
end

# ∂²

function ∂2_2(f::Function, x⃗::Vector{<:Real}, i::Int, h::Float64=1e-3)
  h⃗ = zeros(length(x⃗))
  h⃗[i] = h
  return (f(x⃗ - h⃗) - 2f(x⃗) + f(x⃗ + h⃗)) / h^2
end

function ∂2_4(f::Function, x⃗::Vector{<:Real}, i::Int, h::Float64=1e-3)
  h⃗ = zeros(length(x⃗))
  h⃗[i] = h
  return (-(1 / 12)f(x⃗ - 2h⃗) + (4 / 3)f(x⃗ - h⃗) - (5 / 2)f(x⃗) + (4 / 3)f(x⃗ + h⃗) - (1 / 12)f(x⃗ + 2h⃗)) / h^2
end

function ∂2_6(f::Function, x⃗::Vector{<:Real}, i::Int, h::Float64=1e-3)
  h⃗ = zeros(length(x⃗))
  h⃗[i] = h
  return ((1 / 90)f(x⃗ - 3h⃗) - (3 / 20)f(x⃗ - 2h⃗) + (3 / 2)f(x⃗ - h⃗) - (49 / 18)f(x⃗) + (3 / 2)f(x⃗ + h⃗) - (3 / 20)f(x⃗ + 2h⃗) + (1 / 90)f(x⃗ + 3h⃗)) / h^2
end

function ∂2_8(f::Function, x⃗::Vector{<:Real}, i::Int, h::Float64=1e-3)
  h⃗ = zeros(length(x⃗))
  h⃗[i] = h
  return (-(1 / 560)f(x⃗ - 4h⃗) + (8 / 315)f(x⃗ - 3h⃗) - (1 / 5)f(x⃗ - 2h⃗) + (8 / 5)f(x⃗ - h⃗) - (205 / 72)f(x⃗) + (8 / 5)f(x⃗ + h⃗) - (1 / 5)f(x⃗ + 2h⃗) + (8 / 315)f(x⃗ + 3h⃗) - (1 / 560)f(x⃗ + 4h⃗)) / h^2
end

# ∂³

function ∂3_2(f::Function, x⃗::Vector{<:Real}, i::Int, h::Float64=1e-3)
  h⃗ = zeros(length(x⃗))
  h⃗[i] = h
  return (-(1 / 2)f(x⃗ - 2h⃗) + f(x⃗ - h⃗) - f(x⃗ + h⃗) + (1 / 2)f(x⃗ + 2h⃗)) / h^3
end

function ∂3_4(f::Function, x⃗::Vector{<:Real}, i::Int, h::Float64=1e-3)
  h⃗ = zeros(length(x⃗))
  h⃗[i] = h
  return ((1 / 8)f(x⃗ - 3h⃗) - f(x⃗ - 2h⃗) + (13 / 8)f(x⃗ - h⃗) - (13 / 8)f(x⃗ + h⃗) + f(x⃗ + 2h⃗) - (1 / 8)f(x⃗ + 3h⃗)) / h^3
end

function ∂3_6(f::Function, x⃗::Vector{<:Real}, i::Int, h::Float64=1e-3)
  h⃗ = zeros(length(x⃗))
  h⃗[i] = h
  return (-(7 / 240)f(x⃗ - 4h⃗) + (3 / 10)f(x⃗ - 3h⃗) - (169 / 120)f(x⃗ - 2h⃗) + (61 / 30)f(x⃗ - h⃗) - (61 / 30)f(x⃗ + h⃗) + (169 / 120)f(x⃗ + 2h⃗) - (3 / 10)f(x⃗ + 3h⃗) + (7 / 240)f(x⃗ + 4h⃗)) / h^3
end

# ∂⁴

function ∂4_2(f::Function, x⃗::Vector{<:Real}, i::Int, h::Float64=1e-3)
  h⃗ = zeros(length(x⃗))
  h⃗[i] = h
  return (f(x⃗ - 2h⃗) - 4f(x⃗ - h⃗) + 6f(x⃗) - 4f(x⃗ + h⃗) + f(x⃗ + 2h⃗)) / h^4
end

function ∂4_4(f::Function, x⃗::Vector{<:Real}, i::Int, h::Float64=1e-3)
  h⃗ = zeros(length(x⃗))
  h⃗[i] = h
  return (-(1 / 6)f(x⃗ - 3h⃗) + 2f(x⃗ - 2h⃗) - (13 / 2)f(x⃗ - h⃗) + (28 / 3)f(x⃗) - (13 / 2)f(x⃗ + h⃗) + 2f(x⃗ + 2h⃗) - (1 / 6)f(x⃗ + 3h⃗)) / h^4
end

function ∂4_6(f::Function, x⃗::Vector{<:Real}, i::Int, h::Float64=1e-3)
  h⃗ = zeros(length(x⃗))
  h⃗[i] = h
  return ((7 / 240)f(x⃗ - 4h⃗) - (2 / 5)f(x⃗ - 3h⃗) + (169 / 60)f(x⃗ - 2h⃗) - (122 / 15)f(x⃗ - h⃗) + (91 / 8)f(x⃗) - (122 / 15)f(x⃗ + h⃗) + (169 / 60)f(x⃗ + 2h⃗) - (2 / 5)f(x⃗ + 3h⃗) + (7 / 240)f(x⃗ + 4h⃗)) / h^4
end

# ∂⁵

function ∂5_2(f::Function, x⃗::Vector{<:Real}, i::Int, h::Float64=1e-3)
  h⃗ = zeros(length(x⃗))
  h⃗[i] = h
  return (-(1 / 2)f(x⃗ - 3h⃗) + 2f(x⃗ - 2h⃗) - (5 / 2)f(x⃗ - h⃗) + (5 / 2)f(x⃗ + h⃗) - 2f(x⃗ + 2h⃗) + (1 / 2)f(x⃗ + 3h⃗)) / h^5
end

function ∂5_4(f::Function, x⃗::Vector{<:Real}, i::Int, h::Float64=1e-3)
  h⃗ = zeros(length(x⃗))
  h⃗[i] = h
  return ((1 / 6)f(x⃗ - 4h⃗) - (3 / 2)f(x⃗ - 3h⃗) + (13 / 3)f(x⃗ - 2h⃗) - (29 / 6)f(x⃗ - h⃗) + (29 / 6)f(x⃗ + h⃗) - (13 / 3)f(x⃗ + 2h⃗) + (3 / 2)f(x⃗ + 3h⃗) - (1 / 6)f(x⃗ + 4h⃗)) / h^5
end

function ∂5_6(f::Function, x⃗::Vector{<:Real}, i::Int, h::Float64=1e-3)
  h⃗ = zeros(length(x⃗))
  h⃗[i] = h
  return (-(13 / 288)f(x⃗ - 5h⃗) + (19 / 36)f(x⃗ - 4h⃗) - (87 / 32)f(x⃗ - 3h⃗) + (13 / 2)f(x⃗ - 2h⃗) - (323 / 48)f(x⃗ - h⃗) + (323 / 48)f(x⃗ + h⃗) - (13 / 2)f(x⃗ + 2h⃗) + (87 / 32)f(x⃗ + 3h⃗) - (19 / 36)f(x⃗ + 4h⃗) + (13 / 288)f(x⃗ + 5h⃗)) / h^5
end

# ∂⁶

function ∂6_2(f::Function, x⃗::Vector{<:Real}, i::Int, h::Float64=1e-3)
  h⃗ = zeros(length(x⃗))
  h⃗[i] = h
  return (f(x⃗ - 3h⃗) - 6f(x⃗ - 2h⃗) + 15f(x⃗ - h⃗) - 20f(x⃗) + 15f(x⃗ + h⃗) - 6f(x⃗ + 2h⃗) + f(x⃗ + 3h⃗)) / h^6
end

function ∂6_4(f::Function, x⃗::Vector{<:Real}, i::Int, h::Float64=1e-3)
  h⃗ = zeros(length(x⃗))
  h⃗[i] = h
  return (-(1 / 4)f(x⃗ - 4h⃗) + 3f(x⃗ - 3h⃗) - 13f(x⃗ - 2h⃗) + 29f(x⃗ - h⃗) - (75 / 2)f(x⃗) + 29f(x⃗ + h⃗) - 13f(x⃗ + 2h⃗) + 3f(x⃗ + 3h⃗) - (1 / 4)f(x⃗ + 4h⃗)) / h^6
end

function ∂6_6(f::Function, x⃗::Vector{<:Real}, i::Int, h::Float64=1e-3)
  h⃗ = zeros(length(x⃗))
  h⃗[i] = h
  return ((13 / 240)f(x⃗ - 5h⃗) - (19 / 24)f(x⃗ - 4h⃗) + (87 / 16)f(x⃗ - 3h⃗) - (39 / 2)f(x⃗ - 2h⃗) + (323 / 8)f(x⃗ - h⃗) - (1023 / 20)f(x⃗) + (323 / 8)f(x⃗ - h⃗) - (39 / 2)f(x⃗ - 2h⃗) + (87 / 16)f(x⃗ - 3h⃗) - (19 / 24)f(x⃗ - 4h⃗) + (13 / 240)f(x⃗ - 5h⃗)) / h^6
end
