function ∂1_4_slow(f::Function, x⃗::Vector{<:Real}, i::Int, h::Real=1e-3)::Union{Real,Vector{<:Real}}
  h⃗ = zeros(length(x⃗))
  h⃗[i] = h
  return ((1 / 12)f(x⃗ - 2h⃗) - (2 / 3)f(x⃗ - h⃗) + (2 / 3)f(x⃗ + h⃗) - (1 / 12)f(x⃗ + 2h⃗)) / h
end

function ∂1_4_fast(f::Function, x⃗::Vector{<:Real}, i::Int, h::Real=1e-3)::Union{Real,Vector{<:Real}}
  x_local = copy(x⃗)

  function x(h)
    x_local[i] = x⃗[i]
    x_local[i] += h
    return x_local
  end

  return ((1 / 12)f(x(-2h)) - (2 / 3)f(x(-h)) + (2 / 3)f(x(h)) - (1 / 12)f(x(2h))) / h
end

x⃗ = [Float64(π), Float64(π)]
f(x⃗) = sin(x⃗[1]) + cos(x⃗[2])
