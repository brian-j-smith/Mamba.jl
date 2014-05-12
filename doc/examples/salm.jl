using MCMCsim
using Distributions

## Data
salm = (String => Any)[
  "y" => reshape(
    [15, 21, 29, 16, 18, 21, 16, 26, 33, 27, 41, 60, 33, 38, 41, 20, 27, 42],
     3, 6),
  "x" => [0, 10, 33, 100, 333, 1000],
  "plate" => 3,
  "dose" => 6
]


## Model Specification

model = MCMCModel(

  y = MCMCStochastic(2,
    @modelexpr(alpha, beta, gamma, x, lambda,
      Distribution[
        begin
          mu = exp(alpha + beta * log(x[j] + 10) + gamma * x[j] + lambda[i,j])
          Poisson(mu)
        end
        for i in 1:3, j in 1:6
      ]
    ),
    false
  ),

  alpha = MCMCStochastic(
    :(Normal(0, 1000))
  ),

  beta = MCMCStochastic(
    :(Normal(0, 1000))
  ),

  gamma = MCMCStochastic(
    :(Normal(0, 1000))
  ),

  lambda = MCMCStochastic(2,
    @modelexpr(s2,
      begin
        sigma = sqrt(s2)
        Distribution[Normal(0, sigma) for i in 1:3, j in 1:6]
      end
    ),
    false
  ),

  s2 = MCMCStochastic(
    :(InverseGamma(0.001, 0.001))
  )

)


## Initial Values
inits = [
  ["y" => salm["y"], "alpha" => 0, "beta" => 0, "gamma" => 0, "s2" => 10,
   "lambda" => zeros(3, 6)],
  ["y" => salm["y"], "alpha" => 1, "beta" => 1, "gamma" => 0.01, "s2" => 1,
   "lambda" => zeros(3, 6)]
]


## Sampling Scheme
scheme = [Slice(["alpha", "beta", "gamma"], [1.0, 1.0, 0.1]),
          NUTS(["lambda", "s2"])]
setsamplers!(model, scheme)


## MCMC Simulations
sim = mcmc(model, salm, inits, 10000, burnin=2500, thin=2, chains=2)
describe(sim)
