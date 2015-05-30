using Mamba

include("utils.jl")

test_examples = [
  "blocker",
  "bones",
  "dogs",
  "dyes",
  "epil",
  "equiv",
  "eyes",
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

for t in test_examples
  @everywhere srand(123)
  @runtest "../doc/examples/" t
  gelmandiag(sim) |> show
end
