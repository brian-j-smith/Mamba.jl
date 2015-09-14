include("utils.jl")

test_tutorials = [
  "line"
]

test_samplers = [
  "amm",
  "amwg",
  "bds",
  "nuts",
  "slice"
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

for t in test_extensions
  @everywhere srand(123)
  @runtest "../doc/mcmc/" t
end
