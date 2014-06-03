tutorials = [
  "line"
]

samplers = [
  "amm",
  "amwg",
  "nuts",
  "slice",
  "slicewg"
]

examples = [
  "blocker",
  "bones",
  "dogs",
  "dyes",
  "epil",
  "equiv",
  "inhalers",
#  "kidney",
  "leuk",
  "lsat",
  "magnesium",
#  "mice",
  "oxford",
  "pumps",
  "rats",
  "salm",
  "seeds",
  "stacks",
  "surgical"
]

macro runtest(dir, prefix)
  quote
    fname = $dir * $prefix * ".jl"
    print("\n>>> Testing $fname\n\n")
    include(fname)
  end
end

println("Running tests:")

for t in tutorials
  @runtest "../doc/tutorial/" t
end

for t in samplers
  @runtest "../doc/samplers/" t
end

for t in examples
  @runtest "../doc/examples/" t
end
