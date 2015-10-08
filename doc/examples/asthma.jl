using Mamba

## Data
asthma = Dict{Symbol,Any}(
  :y => 
    [210 60 0 1  1
     88 641 0 4 13
     1    0 0 0  1],
  :M =>
    [272, 746, 2]
)

## Model Specification

model = Model(
  y = Stochastic(2,
    @modelexpr(M,q,
      MultivariateDistribution[
          Multinomial(M[i], vec(q[i,:]))
        for i in 1:length(M)
      ]
    ), false
  ),
  q = Stochastic(2,
    @modelexpr(M,
      MultivariateDistribution[
          Dirichlet(ones(5))
        for i in 1:length(M)
      ]
    ),true
  )
)

## Initial Values
inits = [
  Dict{Symbol,Any}(:y => asthma[:y], :q => vcat([rand(Dirichlet(ones(5)))' for i in 1:3]...))
  for i in 1:3
]

## Sampling Scheme
scheme = [SliceSimplex([:q])]
setsamplers!(model, scheme)


## MCMC Simulations
sim = mcmc(model, asthma, inits, 10000, burnin=2500, thin=1, chains=3)
describe(sim)

