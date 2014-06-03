using MCMCsim
using Distributions

## Model Specification

model = MCMCModel(

  y = MCMCStochastic(1,
    quote
      mu = model[:mu]
      s2 = model[:s2]
      IsoNormal(mu, sqrt(s2))
    end,
    false
  ),

  mu = MCMCLogical(1,
    :(model[:xmat] * model[:beta]),
    false
  ),

  beta = MCMCStochastic(1,
    :(IsoNormal(2, sqrt(1000)))
  ),

  s2 = MCMCStochastic(
    :(InverseGamma(0.001, 0.001))
  )

)

## Hybrid No-U-Turn and Slice Sampling Scheme
scheme1 = [NUTS([:beta]),
           Slice([:s2], [1.0])]

## No-U-Turn Sampling Scheme
scheme2 = [NUTS([:beta, :s2])]

## User-Defined Samplers

Gibbs_beta = MCMCSampler([:beta],
  quote
    beta = model[:beta]
    s2 = model[:s2]
    xmat = model[:xmat]
    y = model[:y]
    beta_mean = mean(beta.distr)
    beta_invcov = invcov(beta.distr)
    Sigma = inv(xmat' * xmat / s2 + beta_invcov)
    mu = Sigma * (xmat' * y / s2 + beta_invcov * beta_mean)
    rand(MvNormal(mu, Sigma))
  end
)

Gibbs_beta = MCMCSampler([:beta],
  @modelexpr(beta, s2, xmat, y,
    begin
      beta_mean = mean(beta.distr)
      beta_invcov = invcov(beta.distr)
      Sigma = inv(xmat' * xmat / s2 + beta_invcov)
      mu = Sigma * (xmat' * y / s2 + beta_invcov * beta_mean)
      rand(MvNormal(mu, Sigma))
    end
  )
)

Gibbs_s2 = MCMCSampler([:s2],
  quote
    beta = model[:beta]
    s2 = model[:s2]
    xmat = model[:xmat]
    y = model[:y]
    a = length(y) / 2.0 + s2.distr.shape
    b = sum((y - xmat * beta).^2) / 2.0 + s2.distr.scale
    rand(InverseGamma(a, b))
  end
)

Gibbs_s2 = MCMCSampler([:s2],
  @modelexpr(beta, s2, xmat, y,
    begin
      a = length(y) / 2.0 + s2.distr.shape
      b = sum((y - xmat * beta).^2) / 2.0 + s2.distr.scale
      rand(InverseGamma(a, b))
    end
  )
)

## User-Defined Sampling Scheme
scheme3 = [Gibbs_beta, Gibbs_s2]

## Sampling Scheme Assignment
setsamplers!(model, scheme1)


## Graph Representation of Nodes
print(graph2dot(model))
graph2dot(model, "lineDAG.dot")


## Data
line = (Symbol => Any)[
  :x => [1, 2, 3, 4, 5],
  :y => [1, 3, 3, 3, 5]
]
line[:xmat] = [ones(5) line[:x]]


## Set Random Number Generator Seed
srand(123)


## Initial Values
inits = [[:y => line[:y],
          :beta => rand(Normal(0, 1), 2),
          :s2 => rand(Gamma(1, 1))]
         for i in 1:3]


## MCMC Simulations

setsamplers!(model, scheme1)
sim1 = mcmc(model, line, inits, 10000, burnin=250, thin=2, chains=3)

setsamplers!(model, scheme2)
sim2 = mcmc(model, line, inits, 10000, burnin=250, thin=2, chains=3)

setsamplers!(model, scheme3)
sim3 = mcmc(model, line, inits, 10000, burnin=250, thin=2, chains=3)


## Brooks, Gelman and Rubin Convergence Diagnostic
gelmandiag(sim1, mpsrf=true, transform=true)

## Geweke Convergence Diagnostic
gewekediag(sim1)


## Summary Statistics
describe(sim1)

## Highest Posterior Density Intervals
hpd(sim1)

## Cross-Correlations
cor(sim1)

## Lag-Autocorrelations
autocor(sim1)

## Deviance Information Criterion
dic(sim1)

## Subset Sampler Output
describe(sim1[1000:5000, ["beta[1]", "beta[2]"], :])


## Development and Testing

setinputs!(model, line)             # Set input node values
setinits!(model, inits[1])          # Set initial values
setsamplers!(model, scheme1)        # Set sampling scheme

showall(model)                      # Show detailed node information

logpdf(model, 1)                    # Log-density sum for block 1
logpdf(model, 2)                    # Block 2
logpdf(model)                       # All blocks

simulate!(model, 1)                 # Simulate draws for block 1
simulate!(model, 2)                 # Block 2
simulate!(model)                    # All blocks
