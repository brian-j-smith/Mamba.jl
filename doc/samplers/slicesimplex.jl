################################################################################
## Multinomial Model
##   y ~ Multinomial(n, rho)
##   rho ~ Dirichlet(1, ..., 1)
################################################################################

using Mamba

## Data
n, k = 100, 5
rho0 = rand(Dirichlet(ones(k)))
y = rand(Multinomial(n, rho0))

## Log-transformed Posterior(rho) + Constant
logf = function(rho::DenseVector)
  logpdf(Multinomial(n, rho), y)
end

## MCMC Simulation with Slice Simplex Sampling
t = 10000
sim = Chains(t, k, names = map(i -> "rho[$i]", 1:k))
rho = SliceSimplexVariate(fill(1 / k, k), logf)
for i in 1:t
  sample!(rho)
  sim[i, :, 1] = rho
end
describe(sim)

p = plot(sim)
draw(p, filename = "slicesimplexplot")
