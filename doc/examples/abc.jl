################################################################################
## Linear Regression
##   y ~ N(b0 + b1 * x, s2)
##   b0, b1 ~ N(0, 1000)
##   s2 ~ invgamma(0.001, 0.001)
################################################################################

using Mamba

## Data
line = Dict{Symbol, Any}(
  :x => [1, 2, 3, 4, 5],
  :y => [1, 3, 3, 3, 5]
)

line[:xmat] = [ones(5) line[:x]]

## Model Specification

model = Model(

  y = Stochastic(1,
    (mu, s2) -> MvNormal(mu, sqrt(s2)),
    false
  ),

  mu = Logical(1,
    (xmat, beta) -> xmat * beta,
    false
  ),

  beta = Stochastic(1,
    () -> MvNormal(2, sqrt(1000))
  ),

  s2 = Stochastic(
    () -> InverseGamma(0.001, 0.001)
  )

)

## Initial Values
inits = [
  Dict{Symbol, Any}(
    :y => line[:y],
    :beta => rand(Normal(0, 1), 2),
    :s2 => rand(Gamma(1, 1))
  )
  for i in 1:3
]

## summary Statistics
summarize = [mean; var; x -> sortperm(x)[1] < sortperm(x)[end]]

epsilon1 = 0.2
sigma1 = 1.0

epsilon2 = 0.3
sigma2 = 8.0

## User-Defined Sampling Scheme
scheme = [ABC([:beta], sigma1, summarize, epsilon1, s = 20, n = 250, kernel = :uniform),
          ABC([:s2],   sigma2, summarize, epsilon2, s = 20, n = 250, kernel = :epanechnikov)]
setsamplers!(model, scheme)

## MCMC Simulation with Approximate Bayesian Computing
sim = mcmc(model, line, inits, 10000, burnin=1000, thin=1, chains=3)
describe(sim)
