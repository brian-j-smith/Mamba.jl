using MCMCsim
using Distributions

## Data
stacks = (String => Any)[
  "y" => [42, 37, 37, 28, 18, 18, 19, 20, 15, 14, 14, 13, 11, 12, 8, 7, 8, 8, 9,
          15, 15],
  "x" =>
    [80 27 89
     80 27 88
     75 25 90
     62 24 87
     62 22 87
     62 23 87
     62 24 93
     62 24 93
     58 23 87
     58 18 80
     58 18 89
     58 17 88
     58 18 82
     58 19 93
     50 18 89
     50 18 86
     50 19 72
     50 19 79
     50 20 80
     56 20 82
     70 20 91]
]
stacks["N"] = size(stacks["x"], 1)
stacks["p"] = size(stacks["x"], 2)

stacks["z"] = hcat(
  ones(21),
  mapslices(x -> (x - mean(x)) / std(x), stacks["x"], 1)
)


## Model Specification

model = MCMCModel(

  y = MCMCStochastic(1,
    @modelexpr(z, beta, s2,
      begin
        mu = z * beta
        sigma = sqrt(s2)
        Distribution[Laplace(mu[i], sigma) for i in 1:length(mu)]
      end
    ),
    false
  ),

  beta = MCMCStochastic(1,
    @modelexpr(p,
      IsoNormal(1 + p, 1000)
    )
  ),

  s2 = MCMCStochastic(
    :(InverseGamma(0.001, 0.001))
  )

)


## Initial Values
inits = [
  ["y" => stacks["y"], "beta" => [10, 0, 0, 0], "s2" => 10],
  ["y" => stacks["y"], "beta" => [1, 1, 1, 1], "s2" => 1]
]


## Sampling Scheme
scheme = [NUTS(["beta"]),
          Slice(["s2"], [1.0])]
setsamplers!(model, scheme)


## MCMC Simulations
sim = mcmc(model, stacks, inits, 10000, burnin=2500, thin=2, chains=2)
describe(sim)
