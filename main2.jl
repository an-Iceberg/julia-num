using CairoMakie

f = Figure()
ax = Axis(
  f[1, 1],
  title="Dodged bars",
  xticks=(1:4, ["left", "middle1", "middle2", "right"])
)

table = (
  cat=[1, 1, 2, 2, 3, 3, 4, 4],
  data=[0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8],
  group=repeat([1, 2], 4),
)
barplot!(
  ax,
  table.cat,
  table.data,
  dodge=table.group,
  color=table.group,
)
save("main2.png", f)
