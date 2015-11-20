################################################################################
## Linear Regression
##   y ~ MvNormal(X * (beta0 .* gamma), 1)
##   gamma ~ DiscreteUniform(0, 1)
################################################################################

using Mamba

## Data
n, p = 25, 10
X = randn(n, p)
beta0 = randn(p)
gamma0 = rand(0:1, p)
y = X * (beta0 .* gamma0) + randn(n)

## Log-transformed Posterior(gamma) + Constant
logf = function(gamma::DenseVector)
  logpdf(MvNormal(X * (beta0 .* gamma), 1.0), y)
end

## MCMC Simulation with Binary MCMC Model Composition
t = 10000
sim = Chains(t, p, names = map(i -> "gamma[$i]", 1:p))
gamma = BMC3Variate(zeros(p))
indexset = collect(combinations(1:p, 1))
for i in 1:t
  bmc3!(gamma, indexset, logf)
  sim[i, :, 1] = gamma
end
describe(sim)

p = plot(sim, [:trace, :mixeddensity])
draw(p, filename = "bmc3plot")
