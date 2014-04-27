using MCMCsim
using Distributions

## Model Specification

line = MCMCModel(

  y = MCMCStochastic(5,
    quote
      mu = model["mu"]
      s2 = model["s2"]
      IsoNormal(mu, sqrt(s2))
    end,
    false
  ),

  mu = MCMCLogical(5,
    :(model["xmat"] * model["beta"]),
    false
  ),

  beta = MCMCStochastic(2,
    :(IsoNormal(2, sqrt(1000)))
  ),

  s2 = MCMCStochastic(
    :(InverseGamma(0.001, 0.001))
  )

)

## Hybrid No-U-Turn and Slice Sampling Scheme
scheme1 = [SamplerNUTS(["beta"]),
           SamplerSlice(["s2"], [1.0])]

## No-U-Turn Sampling Scheme
scheme2 = [SamplerNUTS(["beta", "s2"])]

## User-Defined Samplers

Gibbs_beta = MCMCSampler(["beta"],
  quote
    beta = model["beta"]
    s2 = model["s2"]
    xmat = model["xmat"]
    y = model["y"]
    beta_mean = mean(beta.distr)
    beta_invcov = invcov(beta.distr)
    Sigma = inv(xmat' * xmat / s2 + beta_invcov)
    mu = Sigma * (xmat' * y / s2 + beta_invcov * beta_mean)
    rand(MvNormal(mu, Sigma))
  end
)

Gibbs_s2 = MCMCSampler(["s2"],
  quote
    beta = model["beta"]
    s2 = model["s2"]
    xmat = model["xmat"]
    y = model["y"]
    a = length(y) / 2.0 + s2.distr.shape
    b = sum((y - xmat * beta).^2) / 2.0 + s2.distr.scale
    rand(InverseGamma(a, b))
  end
)

## User-Defined Sampling Scheme
scheme3 = [Gibbs_beta, Gibbs_s2]

## Sampling Scheme Assignment
setsamplers!(line, scheme1)


## Graph Representation of Nodes
print(graph2dot(line))
graph2dot(line, "lineDAG.dot")


## Data
data = (String => Any)[
  "x" => [1, 2, 3, 4, 5],
  "y" => [1, 3, 3, 3, 5]
]
data["xmat"] = [ones(5) data["x"]]


## Initial Values
inits = [["y" => data["y"],
          "beta" => rand(Normal(0, 1), 2),
          "s2" => rand(Gamma(1, 1))]
         for i in 1:3]


## MCMC Simulations

setsamplers!(line, scheme1)
sim1 = mcmc(line, data, inits, 10000, burnin=250, thin=2, chains=3)

setsamplers!(line, scheme2)
sim2 = mcmc(line, data, inits, 10000, burnin=250, thin=2, chains=3)

setsamplers!(line, scheme3)
sim3 = mcmc(line, data, inits, 10000, burnin=250, thin=2, chains=3)


## Brooks, Gelman and Rubin Convergence Diagnostic
gelmandiag(sim1, mpsrf=true, transform=true)


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

setinputs!(line, data)             # Set input node values
setinits!(line, inits[1])          # Set initial values
setsamplers!(line, scheme1)        # Set sampling scheme

showall(line)                      # Show detailed node information

logpdf(line, 1)                    # Log-density sum for block 1
logpdf(line, 2)                    # Block 2
logpdf(line)                       # All blocks

simulate!(line, 1)                 # Simulate draws for block 1
simulate!(line, 2)                 # Block 2
simulate!(line)                    # All blocks
