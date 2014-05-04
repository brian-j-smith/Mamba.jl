using MCMCsim
using Distributions

## Data
pumps = (String => Any)[
  "t" => [94.3, 15.7, 62.9, 126, 5.24, 31.4, 1.05, 1.05, 2.1, 10.5],
  "y" => [5, 1, 5, 14, 3, 19, 1, 1, 4, 22],
  "N" => 10
]


## Model Specification

model = MCMCModel(

  y = MCMCStochastic(1,
    @modelexpr(theta, t, N,
      begin
        lambda = theta .* t
        Distribution[Poisson(lambda[i]) for i in 1:N]
      end
    ),
    false
  ),

  theta = MCMCStochastic(1,
    @modelexpr(alpha, beta, N,
      Distribution[Gamma(alpha, beta) for i in 1:N]
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
scheme = [SamplerSliceWG(["alpha", "beta"], [1.0, 1.0]),
          SamplerAMWG(["theta"], ones(pumps["N"]))]
setsamplers!(model, scheme)


## MCMC Simulations
sim = mcmc(model, pumps, inits, 10000, burnin=2500, thin=2, chains=2)
describe(sim)
