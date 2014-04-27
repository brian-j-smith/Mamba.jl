using MCMCsim
using Distributions

## Data
data = (String => Any)[
  "r" => [10, 23, 23, 26, 17, 5, 53, 55, 32, 46, 10, 8, 10, 8, 23, 0, 3, 22, 15,
          32, 3],
  "n" => [39, 62, 81, 51, 39, 6, 74, 72, 51, 79, 13, 16, 30, 28, 45, 4, 12, 41,
          30, 51, 7],
  "x1" => [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
  "x2" => [0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1],
  "N" => 21
]


## Model Specification

seeds = MCMCModel(

  r = MCMCStochastic(data["N"],
    quote
      eta = model["alpha0"] + model["alpha1"] * model["x1"] +
              model["alpha2"] * model["x2"] +
              model["alpha12"] * model["x1"] .* model["x2"] + model["b"]
      p = 1.0 / (exp(-eta) + 1.0)
      Distribution[Binomial(model["n"][i], p[i]) for i in 1:model["N"]]
    end,
    false
  ),

  b = MCMCStochastic(data["N"],
    :(IsoNormal(model["N"], model["s2"])),
    false
  ),

  alpha0 = MCMCStochastic(
    :(Normal(0.0, 1.0e6))
  ),

  alpha1 = MCMCStochastic(
    :(Normal(0.0, 1.0e6))
  ),

  alpha2 = MCMCStochastic(
    :(Normal(0.0, 1.0e6))
  ),

  alpha12 = MCMCStochastic(
    :(Normal(0.0, 1.0e6))
  ),

  s2 = MCMCStochastic(
    :(InverseGamma(0.001, 0.001))
  )

)


## Initial Values
inits = [
  ["r" => data["r"], "alpha0" => 0, "alpha1" => 0, "alpha2" => 0,
   "alpha12" => 0, "s2" => 0.01, "b" => zeros(data["N"])],
  ["r" => data["r"], "alpha0" => 0, "alpha1" => 0, "alpha2" => 0,
   "alpha12" => 0, "s2" => 1, "b" => zeros(data["N"])]
]


## Sampling Scheme
scheme = [SamplerAMM(["alpha0", "alpha1", "alpha2", "alpha12"], 0.01 * eye(4),
                     adapt=:all),
          SamplerAMWG(["b"], 0.01 * ones(data["N"]), adapt=:all),
          SamplerSlice(["s2"], [1.0])]
setsamplers!(seeds, scheme)


## MCMC Simulations
sim = mcmc(seeds, data, inits, 10000, burnin=2500, thin=2, chains=2)
describe(sim)
