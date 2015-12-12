using Mamba

## Data
pumps = Dict{Symbol, Any}(
  :y => [5, 1, 5, 14, 3, 19, 1, 1, 4, 22],
  :t => [94.3, 15.7, 62.9, 126, 5.24, 31.4, 1.05, 1.05, 2.1, 10.5]
)
pumps[:N] = length(pumps[:y])


## Model Specification
model = Model(

  y = Stochastic(1,
    (theta, t, N) ->
      UnivariateDistribution[
        begin
          lambda = theta[i] * t[i]
          Poisson(lambda)
        end
        for i in 1:N
      ],
    false
  ),

  theta = Stochastic(1,
    (alpha, beta) -> Gamma(alpha, 1 / beta),
    true
  ),

  alpha = Stochastic(
    () -> Exponential(1.0)
  ),

  beta = Stochastic(
    () -> Gamma(0.1, 1.0)
  )

)


## Initial Values
inits = [
  Dict(:y => pumps[:y], :alpha => 1.0, :beta => 1.0,
       :theta => rand(Gamma(1.0, 1.0), pumps[:N])),
  Dict(:y => pumps[:y], :alpha => 10.0, :beta => 10.0,
       :theta => rand(Gamma(10.0, 10.0), pumps[:N]))
]


## Sampling Scheme
scheme = [Slice([:alpha, :beta], 1.0, :univar),
          Slice(:theta, 1.0, :univar)]
setsamplers!(model, scheme)


## MCMC Simulations
sim = mcmc(model, pumps, inits, 10000, burnin=2500, thin=2, chains=2)
describe(sim)


## Posterior Predictive Distribution
ppd = predict(sim, :y)
describe(ppd)
