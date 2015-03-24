include("utils.jl")

tutorials = [
  "line"
]

samplers = [
  "amm",
  "amwg",
  "nuts",
  "slice"
]

println("Running tests:")

for t in tutorials
  @everywhere srand(123)
  @runtest "../doc/tutorial/" t
end

for t in samplers
  @everywhere srand(123)
  @runtest "../doc/samplers/" t
end
