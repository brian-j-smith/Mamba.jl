using Mamba

include("utils.jl")

test_openbugs = [
  "asthma",
  "birats",
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

test_contributed = [
  "line_amwg_slice",
  "pollution"
]

println("Running examples:")

for t in test_openbugs
  @everywhere srand(123)
  @runtest "../doc/examples/" t
  gelmandiag(sim) |> show
end

for t in test_contributed
  @everywhere srand(123)
  @runtest "../doc/examples/" t
end
