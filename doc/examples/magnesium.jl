using MCMCsim
using Distributions

## Data
magnesium = (String => Any)[
  "rt" => [1, 9, 2, 1, 10, 1, 1, 90],
  "nt" => [40, 135, 200, 48, 150, 59, 25, 1159],
  "rc" => [2, 23, 7, 1, 8, 9, 3, 118],
  "nc" => [36, 135, 200, 46, 148, 56, 23, 1157]
]
magnesium["rtx"] = hcat([magnesium["rt"] for i in 1:6]...)'
magnesium["rcx"] = hcat([magnesium["rc"] for i in 1:6]...)'
magnesium["s2"] = 1 / (magnesium["rt"] + 0.5) +
  1 / (magnesium["nt"] - magnesium["rt"] + 0.5) +
  1 / (magnesium["rc"] + 0.5) +
  1 / (magnesium["nc"] - magnesium["rc"] + 0.5)
magnesium["s2_0"] = 1 / mean(1 / magnesium["s2"])


## Model Specification

model = MCMCModel(

  rcx = MCMCStochastic(2,
    @modelexpr(nc, pc,
      Distribution[Binomial(nc[j], pc[i,j]) for i in 1:6, j in 1:8]
    ),
    false
  ),

  pc = MCMCStochastic(2,
    :(Distribution[Uniform(0, 1) for i in 1:6, j in 1:8]),
    false
  ),

  rtx = MCMCStochastic(2,
    @modelexpr(nt, pc, theta,
      Distribution[
        begin
          phi = logit(pc[i,j])
          pt = invlogit(theta[i,j] + phi)
          Binomial(nt[j], pt)
        end
        for i in 1:6, j in 1:8
      ]
    ),
    false
  ),

  theta = MCMCStochastic(2,
    @modelexpr(mu, tau,
      Distribution[Normal(mu[i], tau[i]) for i in 1:6, j in 1:8]
    ),
    false
  ),

  mu = MCMCStochastic(1,
    :(Distribution[Uniform(-10, 10) for i in 1:6]),
    false
  ),

  OR = MCMCLogical(1,
    @modelexpr(mu,
      exp(mu)
    )
  ),

  tau = MCMCLogical(1,
    @modelexpr(priors, s2_0,
      [ sqrt(1 / priors[1]),
        sqrt(priors[2]),
        priors[3],
        sqrt(s2_0 * (1 / priors[4] - 1)),
        sqrt(s2_0) * (1 / priors[5] - 1),
        sqrt(priors[6]) ]
    )
  ),

  priors = MCMCStochastic(1,
    @modelexpr(s2_0,
      Distribution[
        Gamma(0.001, 0.001),
        Uniform(0, 50),
        Uniform(0, 50),
        Uniform(0, 1),
        Uniform(0, 1),
        Truncated(Normal(0, sqrt(s2_0 / erf(0.75))), 0, Inf)
      ]
    ),
    false
  )

)


## Initial Values
inits = [
  ["rcx" => magnesium["rcx"], "rtx" => magnesium["rtx"],
   "theta" => zeros(6, 8), "mu" => fill(-0.5, 6),
   "pc" => fill(0.5, 6, 8), "priors" => [1, 1, 1, 0.5, 0.5, 1]],
  ["rcx" => magnesium["rcx"], "rtx" => magnesium["rtx"],
   "theta" => zeros(6, 8), "mu" => fill(0.5, 6),
   "pc" => fill(0.5, 6, 8), "priors" => [1, 1, 1, 0.5, 0.5, 1]]
]


## Sampling Scheme
scheme = [SamplerAMWG(["theta"], 0.01 * ones(48)),
          SamplerAMM(["mu"], 0.01 * eye(6)),
          SamplerSliceWG(["pc"], fill(0.25, 48)),
          SamplerSliceWG(["priors"], [1.0, 5.0, 5.0, 0.25, 0.25, 5.0])]
setsamplers!(model, scheme)


## MCMC Simulations
sim = mcmc(model, magnesium, inits, 10000, burnin=2500, thin=2, chains=2)
describe(sim)