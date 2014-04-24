################################################################################
## Linear Regression
##   y ~ N(b0 + b1 * x, s2)
##   b0, b1 ~ N(0, 1000)
##   s2 ~ invgamma(0.001, 0.001)
################################################################################

using MCMCsim

## Data
data = [
  "x" => [1, 2, 3, 4, 5],
  "y" => [1, 3, 3, 3, 5]
]

## Log-transformed Posterior(b0, b1, log(s2)) + Constant and Gradient Vector
fx = function(x)
  b0 = x[1]
  b1 = x[2]
  logs2 = x[3]
  r = data["y"] - b0 - b1 * data["x"]
  logf = (-0.5 * length(data["y"]) - 0.001) * logs2 -
           (0.5 * dot(r, r) + 0.001) / exp(logs2) -
           0.5 * b0^2 / 1000 - 0.5 * b1^2 / 1000
  grad = [
    sum(r) / exp(logs2) - b0 / 1000,
    sum(data["x"] .* r) / exp(logs2) - b1 / 1000,
    -0.5 * length(data["y"]) - 0.001 + (0.5 * dot(r, r) + 0.001) / exp(logs2)
  ]
  logf, grad
end

## MCMC Simulation with No-U-Turn Sampling
n = 5000
burnin = 1000
sim = MCMCChains(n - burnin, ["b0", "b1", "s2"], start = (burnin + 1))
theta = VariateNUTS([0.0, 0.0, 0.0])
eps = nutseps(theta, fx)
for i in 1:n
  nuts!(theta, eps, fx, adapt = (i <= burnin))
  if i > burnin
    sim.data[i - burnin,:,1] = [theta[1:2], exp(theta[3])]
  end
end
describe(sim)
