using Mamba, Distributed
@everywhere using Random

include("utils.jl")

const openbugstests = [
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

const contributedtests = [
  "gk",
  "line_abc",
  "line_amwg_slice",
  "pollution"
]

println("Running examples:")

for t in openbugstests
  @everywhere Random.seed!(123)
  @runtest "../doc/examples/" t
  gelmandiag(sim) |> show
end

for t in contributedtests
  @everywhere Random.seed!(123)
  @runtest "../doc/examples/" t
end
