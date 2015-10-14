include("utils.jl")

test_tutorials = [
  "line"
]

test_samplers = [
  "amm",
  "amwg",
  "bmmg",
  "nuts",
  "mala",
  "slice",
  "slicesimplex"
]

test_mcmc = [
  "readcoda"
]

test_extensions = [
  "newunivardist",
  "newmultivardist"
]

println("Running tests:")

for t in test_tutorials
  @everywhere srand(123)
  @runtest "../doc/tutorial/" t
end

for t in test_samplers
  @everywhere srand(123)
  @runtest "../doc/samplers/" t
end

for t in test_mcmc
  @runtest "../doc/mcmc/" t
end

for t in test_extensions
  @everywhere srand(123)
  @runtest "../doc/mcmc/" t
end
