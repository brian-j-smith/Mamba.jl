include("utils.jl")

const tutorialtests = [
  "line"
]

const samplertests = [
  "amm",
  "amwg",
  "bhmc",
  "bmc3",
  "bmg",
  "hmc",
  "mala",
  "nuts",
  "slice",
  "slicesimplex"
]

const mcmctests = [
  "readcoda"
]

const extensiontests = [
  "newunivardist",
  "newmultivardist"
]

println("Running tests:")

for t in tutorialtests
  @everywhere srand(123)
  @runtest "../doc/tutorial/" t
end

for t in samplertests
  @everywhere srand(123)
  @runtest "../doc/samplers/" t
end

for t in mcmctests
  @runtest "../doc/mcmc/" t
end

for t in extensiontests
  @everywhere srand(123)
  @runtest "../doc/mcmc/" t
end
