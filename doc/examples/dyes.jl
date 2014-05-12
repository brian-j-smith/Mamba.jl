using MCMCsim
using Distributions

## Data
dyes = (String => Any)[
  "y" => [
    1545, 1440, 1440, 1520, 1580,
    1540, 1555, 1490, 1560, 1495,
    1595, 1550, 1605, 1510, 1560,
    1445, 1440, 1595, 1465, 1545,
    1595, 1630, 1515, 1635, 1625,
    1520, 1455, 1450, 1480, 1445
  ],
  "batches" => 6,
  "samples" => 5
]

dyes["batch"] = vcat([fill(i, dyes["samples"]) for i in 1:dyes["batches"]]...)
dyes["sample"] = vcat(fill([1:dyes["samples"]], dyes["batches"])...)


## Model Specification

model = MCMCModel(

  y = MCMCStochastic(1,
    @modelexpr(mu, batch, s2_within,
      IsoNormal(mu[batch], sqrt(s2_within))
    ),
    false
  ),

  mu = MCMCStochastic(1,
    @modelexpr(theta, batches, s2_between,
      IsoNormal(theta * ones(batches), sqrt(s2_between))
    ),
    false
  ),

  theta = MCMCStochastic(
    :(Normal(0, 1000))
  ),

  s2_within = MCMCStochastic(
    :(InverseGamma(0.001, 0.001))
  ),

  s2_between = MCMCStochastic(
    :(InverseGamma(0.001, 0.001))
  )

)


## Initial Values
inits = [
  ["y" => dyes["y"], "theta" => 1500, "s2_within" => 1, "s2_between" => 1,
   "mu" => fill(1500, dyes["batches"])],
  ["y" => dyes["y"], "theta" => 3000, "s2_within" => 10, "s2_between" => 10,
   "mu" => fill(3000, dyes["batches"])]
]


## Sampling Scheme
scheme = [NUTS(["mu", "theta"]),
          SliceWG(["s2_within", "s2_between"], [1000.0, 1000.0])]
setsamplers!(model, scheme)


## MCMC Simulations
sim = mcmc(model, dyes, inits, 10000, burnin=2500, thin=2, chains=2)
describe(sim)
