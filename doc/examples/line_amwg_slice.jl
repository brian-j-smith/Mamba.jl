################################################################################
## Linear Regression
##   y ~ N(b0 + b1 * x, s2)
##   b0, b1 ~ N(0, 1000)
##   s2 ~ invgamma(0.001, 0.001)
################################################################################

using Mamba

## Data
data = Dict{Symbol, Any}(
  :x => [1, 2, 3, 4, 5],
  :y => [1, 3, 3, 3, 5]
)

## Log-transformed unnormalized joint posterior for b0, b1, and log(s2)
logf = function(x::DenseVector)
   b0 = x[1]
   b1 = x[2]
   logs2 = x[3]
   (-0.5 * length(data[:y]) - 0.001) * logs2 -
     (0.5 * sum(abs2, data[:y] .- b0 .- b1 .* data[:x]) + 0.001) / exp(logs2) -
     0.5 * b0^2 / 1000 - 0.5 * b1^2 / 1000
end

## Log-transformed unnormalized full conditional densities for the model
## parameters beta and log(s2) defined below in the MCMC simulation
logf_beta(x) = logf([x; logs2])
logf_logs2(x) = logf([beta; x])

## MCMC simulation
n = 10000
burnin = 1000
sim = Chains(n, 3, names = ["b0", "b1", "s2"])
beta = AMWGVariate([0.0, 0.0], 1.0, logf_beta)
logs2 = SliceMultivariate([0.0], 5.0, logf_logs2)
for i in 1:n
  sample!(beta, adapt = (i <= burnin))
  sample!(logs2)
  sim[i, :, 1] = [beta; exp.(logs2)]
end
describe(sim)
