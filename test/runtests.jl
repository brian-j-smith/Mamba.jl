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
  @runtest "../doc/tutorial/" t
end

for t in samplers
  @runtest "../doc/samplers/" t
end
