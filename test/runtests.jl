using Distributed
using Random
using LinearAlgebra

include("utils.jl")

const tutorialtests = [
  "line"
]

const samplertests = [
  "amm",
  "amwg",
  "bhmc",
  "bia",
  "bmc3",
  "bmg",
  "hmc",
  "mala",
  "nuts",
  "rwm",
  "slice",
  "slicesimplex"
]

const mcmctests = [
  "discretediag",
  "readcoda"
]

const extensiontests = [
  "newunivardist",
  "newmultivardist"
]

println("Running tests:")

for t in tutorialtests
  @everywhere Random.seed!(123)
  @runtest "../doc/tutorial/" t
end

for t in samplertests
  @everywhere Random.seed!(123)
  @runtest "../doc/samplers/" t
end

for t in mcmctests
  @runtest "../doc/mcmc/" t
end

for t in extensiontests
  @everywhere Random.seed!(123)
  @runtest "../doc/mcmc/" t
end
