## Define a new univariate Distribution type for Mamba.
## The definition must be placed within an unevaluated quote block.
@everywhere extensions = quote

  ## Load needed packages and import methods to be extended
  using Distributions
  import Distributions: minimum, maximum, logpdf

  ## Type declaration
  type NewUnivarDist <: ContinuousUnivariateDistribution
    mu::Float64
    sigma::Float64
  end

  ## The following method functions must be implemented

  ## Minimum and maximum support values
  minimum(d::NewUnivarDist) = -Inf
  maximum(d::NewUnivarDist) = Inf

  ## Normalized or unnormalized log-density value
  function logpdf(d::NewUnivarDist, x::Real)
    -log(d.sigma) - 0.5 * ((x - d.mu) / d.sigma)^2
  end

end

## Test the extensions in a temporary module (optional)
module Testing end
eval(Testing, extensions)
d = Testing.NewUnivarDist(0.0, 1.0)
Testing.minimum(d)
Testing.maximum(d)
Testing.insupport(d, 2.0)
Testing.logpdf(d, 2.0)

## Add the extensions
using Mamba
@everywhere eval(extensions)

## Implement a Mamba model using the new distribution
model = Model(

  y = Stochastic(1,
    (mu, s2) ->
      begin
        sigma = sqrt(s2)
        UnivariateDistribution[
          NewUnivarDist(mu[i], sigma) for i in 1:length(mu)
        ]
      end,
    false
  ),

  mu = Logical(1,
    (xmat, beta) ->
      xmat * beta,
    false
  ),

  beta = Stochastic(1,
    () -> MvNormal(2, sqrt(1000))
  ),

  s2 = Stochastic(
    () -> InverseGamma(0.001, 0.001)
  )

)

## Sampling Scheme
scheme = [NUTS(:beta),
          Slice(:s2, 3.0)]

## Sampling Scheme Assignment
setsamplers!(model, scheme)

## Data
line = Dict{Symbol, Any}(
  :x => [1, 2, 3, 4, 5],
  :y => [1, 3, 3, 3, 5]
)
line[:xmat] = [ones(5) line[:x]]

## Initial Values
inits = [
  Dict{Symbol, Any}(
    :y => line[:y],
    :beta => rand(Normal(0, 1), 2),
    :s2 => rand(Gamma(1, 1))
  )
  for i in 1:3
]

## MCMC Simulation
sim = mcmc(model, line, inits, 10000, burnin=250, thin=2, chains=3)
describe(sim)
