using MCMCsim
using Distributions

## Data
pumps = (String => Any)[
  "y" => [5, 1, 5, 14, 3, 19, 1, 1, 4, 22],
  "t" => [94.3, 15.7, 62.9, 126, 5.24, 31.4, 1.05, 1.05, 2.1, 10.5]
]
pumps["N"] = length(pumps["y"])


## Model Specification

model = MCMCModel(

  y = MCMCStochastic(1,
    @modelexpr(theta, t, N,
      Distribution[
        begin
          lambda = theta[i] * t[i]
          Poisson(lambda)
        end
        for i in 1:N
      ]
    ),
    false
  ),

  theta = MCMCStochastic(1,
    @modelexpr(alpha, beta,
      Gamma(alpha, 1 / beta)
    )
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
  ["y" => pumps["y"], "alpha" => 1.0, "beta" => 1.0,
   "theta" => rand(Gamma(1.0, 1.0), pumps["N"])],
  ["y" => pumps["y"], "alpha" => 10.0, "beta" => 10.0,
   "theta" => rand(Gamma(10.0, 10.0), pumps["N"])]
]


## Sampling Scheme
scheme = [SliceWG(["alpha", "beta"], [1.0, 1.0]),
          AMWG(["theta"], ones(pumps["N"]))]
setsamplers!(model, scheme)


## MCMC Simulations
sim = mcmc(model, pumps, inits, 10000, burnin=2500, thin=2, chains=2)
describe(sim)
