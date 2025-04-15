using Format

nums = [1000, 20_000, 300_00, 4_000_000, 50_000_000]

nums_strings = [replace(cfmt("%'d", n), "," => "'") for n in nums]
# nums_strings = [replace(n_s, "," => "'") for n_s in nums_strings]
# nums_strings = [join([format("{} ", c) for c in s]) for s in nums_strings]

display(nums_strings)

println([format("{:.0e}", n) for n in nums])

nums_strings = []

println(nums)

for n in nums
  str = ""

  if n >= 1_000 && n < 1_000_000
    str = format("{:d}k", n / 1_000)
  elseif n >= 1_000_000 && n < 1_000_000_000
    str = format("{:d}M", n / 1_000_000)
  end

  push!(nums_strings, str)
end

println(nums_strings)
