################################################################################
## Linear Regression
##   y ~ N(b0 + b1 * x, s2)
##   b0, b1 ~ N(0, 1000)
##   s2 ~ invgamma(0.001, 0.001)
################################################################################

using MCMCsim

## Data
data = [
  :x => [1, 2, 3, 4, 5],
  :y => [1, 3, 3, 3, 5]
]

## Log-transformed Posterior(b0, b1, log(s2)) + Constant
logf = function(x)
   b0 = x[1]
   b1 = x[2]
   logs2 = x[3]
   r = data[:y] .- b0 .- b1 * data[:x]
   (-0.5 * length(data[:y]) - 0.001) * logs2 -
     (0.5 * dot(r, r) + 0.001) / exp(logs2) -
     0.5 * b0^2 / 1000 - 0.5 * b1^2 / 1000
end

## MCMC Simulation with Adaptive Metopolis-within-Gibbs Sampling
n = 5000
burnin = 1000
sim = MCMCChains(n, ["b0", "b1", "s2"])
theta = VariateAMWG([0.0, 0.0, 0.0])
sigma = ones(3)
for i in 1:n
  amwg!(theta, sigma, logf, adapt = (i <= burnin))
  sim.value[i,:,1] = [theta[1:2], exp(theta[3])]
end
describe(sim)
