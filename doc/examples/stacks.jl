using Distributed
@everywhere using Mamba, LinearAlgebra

## Data
stacks = Dict{Symbol, Any}(
  :y => [42, 37, 37, 28, 18, 18, 19, 20, 15, 14, 14, 13, 11, 12, 8, 7, 8, 8, 9,
         15, 15],
  :x =>
    [80 27 89
     80 27 88
     75 25 90
     62 24 87
     62 22 87
     62 23 87
     62 24 93
     62 24 93
     58 23 87
     58 18 80
     58 18 89
     58 17 88
     58 18 82
     58 19 93
     50 18 89
     50 18 86
     50 19 72
     50 19 79
     50 20 80
     56 20 82
     70 20 91]
)
stacks[:N] = size(stacks[:x], 1)
stacks[:p] = size(stacks[:x], 2)

stacks[:meanx] = map(j -> mean(stacks[:x][:, j]), 1:stacks[:p])
stacks[:sdx] = map(j -> std(stacks[:x][:, j]), 1:stacks[:p])
stacks[:z] = Float64[
  (stacks[:x][i, j] - stacks[:meanx][j]) / stacks[:sdx][j]
  for i in 1:stacks[:N], j in 1:stacks[:p]
]


## Model Specification
model = Model(

  y = Stochastic(1,
    (mu, s2, N) ->
      UnivariateDistribution[Laplace(mu[i], s2) for i in 1:N],
    false
  ),

  beta0 = Stochastic(
    () -> Normal(0, 1000),
    false
  ),

  beta = Stochastic(1,
    () -> Normal(0, 1000),
    false
  ),

  mu = Logical(1,
    (beta0, z, beta) -> beta0 .+ z * beta,
    false
  ),

  s2 = Stochastic(
    () -> InverseGamma(0.001, 0.001),
    false
  ),

  sigma = Logical(
    s2 -> sqrt(2.0) * s2
  ),

  b0 = Logical(
    (beta0, b, meanx) -> beta0 - dot(b, meanx)
  ),

  b = Logical(1,
    (beta, sdx) -> beta ./ sdx
  ),

  outlier = Logical(1,
    (y, mu, sigma, N) ->
      Float64[abs((y[i] - mu[i]) / sigma) > 2.5 for i in 1:N],
    [1, 3, 4, 21]
  )

)


## Initial Values
inits = [
  Dict(:y => stacks[:y], :beta0 => 10, :beta => [0, 0, 0], :s2 => 10),
  Dict(:y => stacks[:y], :beta0 => 1, :beta => [1, 1, 1], :s2 => 1)
]


## Sampling Scheme
scheme = [NUTS([:beta0, :beta]),
          Slice(:s2, 1.0)]
setsamplers!(model, scheme)


## MCMC Simulations
sim = mcmc(model, stacks, inits, 10000, burnin=2500, thin=2, chains=2)
describe(sim)
