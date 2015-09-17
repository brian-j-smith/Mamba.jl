################################################################################
## Multinomial Model
##   y ~ Multinomial(n,rho)
##   rho ~ Dirichlet(1,...,1)
################################################################################

using Mamba

k = 5
n = 200
rho0 = rand(Dirichlet(ones(k)))

y = rand(Multinomial(n,rho0))

## Log-transformed Posterior(rho) + Constant
logf = function(rho)
  logpdf(Multinomial(n,rho), y)
end


## MCMC Simulation with Slice Simplex Sampling
t = 10000
sim = Chains(t, k, names = map(i -> "rho[$i]", 1:k))
rho = SliceSimplexVariate(ones(k)/k)
for i in 1:t
  slicesimplex!(rho, 1.0, logf)
  sim[i,:,1] = rho
end
describe(sim)

p = plot(sim)
draw(p, filename = "slicesimplexplot")