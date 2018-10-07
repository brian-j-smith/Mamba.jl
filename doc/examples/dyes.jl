using Mamba

## Data
dyes = Dict{Symbol, Any}(
  :y =>
    [1545, 1440, 1440, 1520, 1580,
     1540, 1555, 1490, 1560, 1495,
     1595, 1550, 1605, 1510, 1560,
     1445, 1440, 1595, 1465, 1545,
     1595, 1630, 1515, 1635, 1625,
     1520, 1455, 1450, 1480, 1445],
  :batches => 6,
  :samples => 5
)

dyes[:batch] = vcat([fill(i, dyes[:samples]) for i in 1:dyes[:batches]]...)
dyes[:sample] = vcat(fill(collect(1:dyes[:samples]), dyes[:batches])...)


## Model Specification

model = Model(

  y = Stochastic(1,
    (mu, batch, s2_within) -> MvNormal(mu[batch], sqrt(s2_within)),
    false
  ),

  mu = Stochastic(1,
    (theta, batches, s2_between) -> Normal(theta, sqrt(s2_between))
  ),

  theta = Stochastic(
    () -> Normal(0, 1000)
  ),

  s2_within = Stochastic(
    () -> InverseGamma(0.001, 0.001)
  ),

  s2_between = Stochastic(
    () -> InverseGamma(0.001, 0.001)
  )

)


## Initial Values
inits = [
  Dict(:y => dyes[:y], :theta => 1500, :s2_within => 1, :s2_between => 1,
       :mu => fill(1500, dyes[:batches])),
  Dict(:y => dyes[:y], :theta => 3000, :s2_within => 10, :s2_between => 10,
       :mu => fill(3000, dyes[:batches]))
]


## Sampling Schemes
scheme = [NUTS([:mu, :theta]),
          Slice([:s2_within, :s2_between], 1000.0)]

scheme2 = [MALA(:theta, 50.0),
           MALA(:mu, 50.0, Matrix{Float64}(I, dyes[:batches], dyes[:batches])),
           Slice([:s2_within, :s2_between], 1000.0)]

scheme3 = [HMC(:theta, 10.0, 5),
           HMC(:mu, 10.0, 5, Matrix{Float64}(I, dyes[:batches], dyes[:batches])),
           Slice([:s2_within, :s2_between], 1000.0)]

scheme4 = [RWM(:theta, 50.0, proposal=Cosine),
           RWM(:mu, 50.0),
           Slice([:s2_within, :s2_between], 1000.0)]


## MCMC Simulations
setsamplers!(model, scheme)
sim = mcmc(model, dyes, inits, 10000, burnin=2500, thin=2, chains=2)
describe(sim)

setsamplers!(model, scheme2)
sim2 = mcmc(model, dyes, inits, 10000, burnin=2500, thin=2, chains=2)
describe(sim2)

setsamplers!(model, scheme3)
sim3 = mcmc(model, dyes, inits, 10000, burnin=2500, thin=2, chains=2)
describe(sim3)

setsamplers!(model, scheme4)
sim4 = mcmc(model, dyes, inits, 10000, burnin=2500, thin=2, chains=2)
describe(sim4)
