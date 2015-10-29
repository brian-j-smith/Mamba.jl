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
summarize = [mean; var; findminat; x -> convert(Int64, sortperm(x)[1] < sortperm(x)[end])]

## Distance function
rho = (x, y) -> sqrt(sum((x .- y).^2))/sqrt(sum(x.^2))

epsilon = 0.25
sigma = 5.0

## User-Defined Sampling Scheme
scheme = [ABC([:beta, :s2], sigma, summarize, epsilon, rho = rho, s = 20, n = 250, kernel = :uniform)]
setsamplers!(model, scheme)

## MCMC Simulation with Approximate Bayesian Computing
sim = mcmc(model, line, inits, 10000, burnin=1000, thin=1, chains=3)
describe(sim)
