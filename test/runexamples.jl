include("utils.jl")

examples = [
  "blocker",
  "bones",
  "dogs",
  "dyes",
  "epil",
  "equiv",
  "inhalers",
  "leuk",
  "lsat",
  "magnesium",
  "mice",
  "oxford",
  "pumps",
  "rats",
  "salm",
  "seeds",
  "stacks",
  "surgical"
]

println("Running examples:")

srand(123)

for t in examples
  @runtest "../doc/examples/" t
end
