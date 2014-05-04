using MCMCsim
using Distributions

## Data
equiv = (String => Any)[
  "group" => [1, 1, 2, 2, 2, 1, 1, 1, 2, 2],
  "y" => reshape(
    [1.40, 1.65,
     1.64, 1.57,
     1.44, 1.58,
     1.36, 1.68,
     1.65, 1.69,
     1.08, 1.31,
     1.09, 1.43,
     1.25, 1.44,
     1.25, 1.39,
     1.30, 1.52], 2, 10)'
]
equiv["T"] = [equiv["group"] 3 - equiv["group"]]


## Model Specification

model = MCMCModel(

  y = MCMCStochastic(2,
    @modelexpr(delta, mu, phi, pi, s2_1, T,
      begin
        sigma = sqrt(s2_1)
        Distribution[
          begin
            m = mu + (-1)^(T[i,j]-1) * phi / 2 + (-1)^(j-1) * pi / 2 +
                delta[i,j]
            Normal(m, sigma)
          end
          for i in 1:10, j in 1:2
        ]
      end
    ),
    false
  ),

  delta = MCMCStochastic(2,
    @modelexpr(s2_2,
      begin
        sigma = sqrt(s2_2)
        Distribution[Normal(0, sigma) for i in 1:10, j in 1:2]
      end
    ),
    false
  ),

  mu = MCMCStochastic(
    :(Normal(0, 1e6))
  ),

  phi = MCMCStochastic(
    :(Normal(0, 1e6))
  ),

  theta = MCMCLogical(
    @modelexpr(phi,
      exp(phi)
    )
  ),

  pi = MCMCStochastic(
    :(Normal(0, 1e6))
  ),

  s2_1 = MCMCStochastic(
    :(InverseGamma(0.001, 0.001))
  ),

  s2_2 = MCMCStochastic(
    :(InverseGamma(0.001, 0.001))
  ),

  model = MCMCLogical(
    @modelexpr(theta,
      int(0.8 < theta < 1.2)
    )
  )

)


## Initial Values
inits = [
  ["y" => equiv["y"], "delta" => zeros(10,2), "mu" => 0, "phi" => 0,
   "pi" => 0, "s2_1" => 1, "s2_2" => 1],
  ["y" => equiv["y"], "delta" => zeros(10,2), "mu" => 10, "phi" => 10,
   "pi" => 10, "s2_1" => 10, "s2_2" => 10]
]


## Sampling Scheme
scheme = [SamplerNUTS(["delta"]),
          SamplerAMWG(["mu", "phi", "pi"], 0.1 * ones(3)),
          SamplerSliceWG(["s2_1", "s2_2"], ones(2))]
setsamplers!(model, scheme)


## MCMC Simulations
sim = mcmc(model, equiv, inits, 10000, burnin=2500, thin=2, chains=2)
describe(sim)
