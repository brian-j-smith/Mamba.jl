using Mamba
using Distributions

## Data
surgical = (Symbol => Any)[
  :r => [0, 18, 8, 46, 8, 13, 9, 31, 14, 8, 29, 24],
  :n => [47, 148, 119, 810, 211, 196, 148, 215, 207, 97, 256, 360]
]
surgical[:N] = length(surgical[:r])


## Model Specification

model = Model(

  r = Stochastic(1,
    @modelexpr(n, p, N,
      Distribution[Binomial(n[i], p[i]) for i in 1:N]
    ),
    false
  ),

  p = Logical(1,
    @modelexpr(b,
      invlogit(b)
    )
  ),

  b = Stochastic(1,
    @modelexpr(mu, s2,
      Normal(mu, sqrt(s2))
    ),
    false
  ),

  mu = Stochastic(
    :(Normal(0, 1000))
  ),

  pop_mean = Logical(
    @modelexpr(mu,
      invlogit(mu)
    )
  ),

  s2 = Stochastic(
    :(InverseGamma(0.001, 0.001))
  )

)

srand(123)

## Initial Values
inits = [
  [:r => surgical[:r], :b => fill(0.1, surgical[:N]), :s2 => 1, :mu => 0],
  [:r => surgical[:r], :b => fill(0.5, surgical[:N]), :s2 => 10, :mu => 1]
]


## Sampling Scheme
scheme = [NUTS([:b]),
          Slice([:mu, :s2], [1.0, 1.0])]
setsamplers!(model, scheme)


## MCMC Simulations
sim = mcmc(model, surgical, inits, 10000, burnin=2500, thin=2, chains=2)
describe(sim)
