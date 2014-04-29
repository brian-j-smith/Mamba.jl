using MCMCsim
using Distributions

## Data
data = (String => Any)[
  "r" => [0, 18, 8, 46, 8, 13, 9, 31, 14, 8, 29, 24],
  "n" => [47, 148, 119, 810, 211, 196, 148, 215, 207, 97, 256, 360],
  "N" => 12
]


## Model Specification

surgical = MCMCModel(

  r = MCMCStochastic(data["N"],
    quote
      Distribution[Binomial(model["n"][i], model["p"][i]) for i in 1:model["N"]]
    end,
    false
  ),

  b = MCMCStochastic(data["N"],
    :(IsoNormal(model["mu"] * ones(model["N"]), model["s2"])),
    false
  ),

  p = MCMCLogical(data["N"],
    :(1.0 / (exp(-model["b"]) + 1.0))
  ),

  mu = MCMCStochastic(
    :(Normal(0.0, 1.0e6))
  ),

  pop_mean = MCMCLogical(
    :(1.0 / (exp(-model["mu"]) + 1.0))
  ),

  s2 = MCMCStochastic(
    :(InverseGamma(0.001, 0.001))
  )

)


## Initial Values
inits = [
  ["r" => data["r"], "b" => fill(0.1, data["N"]), "s2" => 1, "mu" => 0],
  ["r" => data["r"], "b" => fill(0.5, data["N"]), "s2" => 10, "mu" => 1]
]


## Sampling Scheme
scheme = [SamplerNUTS(["b"]),
          SamplerSlice(["mu", "s2"], [1.0, 1.0])]
setsamplers!(surgical, scheme)


## MCMC Simulations
sim = mcmc(surgical, data, inits, 10000, burnin=2500, thin=2, chains=2)
describe(sim)
