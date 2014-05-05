using MCMCsim
using Distributions

## Data
mice = (String => Any)[
  "t" => reshape(
    [12, 1, 21, 25, 11, 26, 27, 30, 13, 12, 21, 20, 23, 25, 23, 29, 35, NaN, 31,
     36, 32, 27, 23, 12, 18, NaN, NaN, 38, 29, 30, NaN, 32, NaN, NaN, NaN, NaN,
     25, 30, 37, 27, 22, 26, NaN, 28, 19, 15, 12, 35, 35, 10, 22, 18, NaN, 12,
     NaN, NaN, 31, 24, 37, 29, 27, 18, 22, 13, 18, 29, 28, NaN, 16, 22, 26, 19,
     NaN, NaN, 17, 28, 26, 12, 17, 26], 20, 4)',
  "tcensor" => reshape(
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 40, 0, 0, 0, 0, 0, 0, 0,
     40, 40, 0, 0, 0, 40, 0, 40, 40, 40, 40, 0, 0, 0, 0, 0, 0, 10, 0, 0, 0, 0,
     0, 0, 0, 0, 0, 24, 0, 40, 40, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 20, 0, 0, 0,
     0, 29, 10, 0, 0, 0, 0, 0, 0], 20, 4)',
   "M" => 4,
   "N" => 20
]


## Model Specification

model = MCMCModel(

  t = MCMCStochastic(2,
    @modelexpr(tau, beta, tcensor, M, N,
      Distribution[
        Truncated(Weibull(tau, exp(-beta[i])), tcensor[i,j], Inf)
        for i in 1:M, j in 1:N
      ]
    ),
    false
  ),

  tau = MCMCStochastic(
    :(Gamma(1, 0.001))
  ),

  beta = MCMCStochastic(1,
    @modelexpr(M,
      IsoNormal(M, 100)
    )
  )

)


## Initial Values
inits = [
  ["t" => mice["t"], "tcensor" => mice["tcensor"],
   "beta" => -ones(mice["M"]), "tau" => 1],
  ["t" => mice["t"], "tcensor" => mice["tcensor"],
   "beta" => zeros(mice["M"]), "tau" => 2]
]


## Sampling Scheme
scheme = []
setsamplers!(model, scheme)


## MCMC Simulations
sim = mcmc(model, mice, inits, 10000, burnin=2500, thin=2, chains=2)
describe(sim)
