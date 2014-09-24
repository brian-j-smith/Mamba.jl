using Mamba

include("utils.jl")

examples = [
  "blocker",
  "bones",
  "dogs",
  "dyes",
  "epil",
  "equiv",
  "inhalers",
  "jaws",
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

for t in examples
  @everywhere srand(123)
  @runtest "../doc/examples/" t
  gelmandiag(sim) |> show
end
