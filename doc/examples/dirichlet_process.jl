#!/Applications/Julia-0.6.app/Contents/Resources/julia/bin/julia



## Define a new univariate Distribution type for Mamba.
## The definition must be placed within an unevaluated quote block.

## Data
data = Dict{Symbol, Any}(
  :y => randn(100) + 5
)

## Maximum number of active components
number_of_components = 40


function stick_breaking(beta)
  cum_prod = 1 - beta[1]
  weights = ones(number_of_components)
  weights[1] = beta[1]
  for i = 2:number_of_components
    weights[i] = beta[i]*cum_prod
    cum_prod = cum_prod*(1-beta[i])
  end
  weights
end




@everywhere extensions = quote

  ## Load needed packages and import methods to be extended
  using Distributions
  import Distributions: minimum, maximum, logpdf

  ## Type declaration
  type NewUnivarDist <: ContinuousUnivariateDistribution
    mu::Array{Float64,1}
    sigma::Float64
    w::Array{Float64,1}
  end

  ## The following method functions must be implemented

  ## Minimum and maximum support values
  minimum(d::NewUnivarDist) = -Inf
  maximum(d::NewUnivarDist) = Inf

  ## Normalized or unnormalized log-density value
  function logpdf(d::NewUnivarDist, x::Real)
    logdens = 0
    for i in 1:40
      logdens = logdens  +  d.w[i]*(1.0/sqrt((2*pi*d.sigma)))*exp( - 0.5 * ((x - d.mu[i]) / d.sigma)^2)
    end
    log(logdens)
  end


end


## Add the extensions
using Mamba
@everywhere eval(extensions)
## Implement a Mamba model using the new distribution
model = Model(

  y = Stochastic(1,
    (beta,mu,tau) ->
      begin
        weights = stick_breaking(beta)
        UnivariateDistribution[
          NewUnivarDist(mu, tau,weights) for i in 1:100
        ]
      end,
    false
  ),
  ytilde = Logical(1,
    (beta,mu,tau) ->
      begin
        weights = stick_breaking(beta)
        idx = wsample(1:number_of_components,weights )
        rands = rand(Normal(mu[idx], tau)
        ,1)
        rands
      end
  ),
  beta = Stochastic(1,
    (alpha) ->
     Beta(1, alpha)
    ),

  mu = Stochastic(1,
    (tau) ->
    begin
      Normal(0,tau)
    end
  ),

  tau = Stochastic(
    () -> InverseGamma(1,1)
  ),
  alpha = Stochastic(
    () -> Uniform(0,1)
  )
)


## Sampling Scheme
## Initial Values
inits = [
  Dict{Symbol, Any}(
    :y => data[:y],
    :beta =>rand(number_of_components)*1,
    :mu => rand(number_of_components),
    :tau => 1,
    :alpha => .1,
  ) for i in 1:3
]


## Sampling Scheme

scheme = [Slice([:beta,:mu,:tau,:alpha],3.0)]

setsamplers!(model, scheme)


## MCMC Simulations
sim = mcmc(model, data, inits, 10000, burnin=2500, thin=2, chains=3)
describe(sim)



