using MCMCsim
using Distributions

## Data
data = (String => Any)[
  "t" => [94.3, 15.7, 62.9, 126, 5.24, 31.4, 1.05, 1.05, 2.1, 10.5],
  "y" => [5, 1, 5, 14, 3, 19, 1, 1, 4, 22]
]
data["N"] = length(data["t"])


## Model Specification

pumps = MCMCModel(

  theta = MCMCStochastic(data["N"],
    quote
      alpha = model["alpha"]
      beta = model["beta"]
      Distribution[Gamma(alpha, beta) for i in 1:model["N"]]
    end
  ),

  y = MCMCStochastic(data["N"],
    quote
      lambda = model["theta"] .* model["t"]
      Distribution[Poisson(lambda[i]) for i in 1:model["N"]]
    end,
    false
  ),

  alpha = MCMCStochastic(
    :(Exponential(1.0))
  ),

  beta = MCMCStochastic(
    :(Gamma(0.1, 1.0))
  )

)


## Initial Values
inits = [
  ["y" => data["y"], "alpha" => 1.0, "beta" => 1.0,
   "theta" => rand(Gamma(1.0, 1.0), data["N"])],
  ["y" => data["y"], "alpha" => 10.0, "beta" => 10.0,
   "theta" => rand(Gamma(10.0, 10.0), data["N"])]
]


## Sampling Scheme
scheme = [SamplerSliceWG(["alpha", "beta"], [1.0, 1.0]),
          SamplerAMWG(["theta"], ones(data["N"]))]
setsamplers!(pumps, scheme)


## MCMC Simulations
sim = mcmc(pumps, data, inits, 10000, burnin=2500, thin=2, chains=2)
describe(sim)
