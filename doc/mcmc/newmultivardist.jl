## Define a new multivariate Distribution type for Mamba.
## The definition must be placed within an unevaluated quote block.
@everywhere extensions = quote

  ## Load needed packages and import methods to be extended
  using Distributions
  import Distributions: length, insupport, _logpdf

  ## Type declaration
  mutable struct NewMultivarDist <: ContinuousMultivariateDistribution
    mu::Vector{Float64}
    sigma::Float64
  end

  ## The following method functions must be implemented

  ## Dimension of the distribution
  length(d::NewMultivarDist) = length(d.mu)

  ## Logical indicating whether x is in the support
  function insupport(d::NewMultivarDist, x::AbstractVector{T}) where {T<:Real}
    length(d) == length(x) && all(isfinite.(x))
  end

  ## Normalized or unnormalized log-density value
  function _logpdf(d::NewMultivarDist, x::AbstractVector{T}) where {T<:Real}
    -length(x) * log(d.sigma) - 0.5 * sum(abs2, x - d.mu) / d.sigma^2
  end

end

## Test the extensions in a temporary module (optional)
module Testing end
Core.eval(Testing, extensions)
d = Testing.NewMultivarDist([0.0, 0.0], 1.0)
Testing.insupport(d, [2.0, 3.0])
Testing.logpdf(d, [2.0, 3.0])

## Add the extensions
using Mamba
@everywhere eval(extensions)

## Implement a Mamba model using the new distribution
model = Model(

  y = Stochastic(1,
    (mu, s2) -> NewMultivarDist(mu, sqrt(s2)),
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
