using Distributed
@everywhere using Mamba

## Data
lsat = Dict{Symbol, Any}(
  :culm =>
    [3, 9, 11, 22, 23, 24, 27, 31, 32, 40, 40, 56, 56, 59, 61, 76, 86, 115, 129,
     210, 213, 241, 256, 336, 352, 408, 429, 602, 613, 674, 702, 1000],
  :response =>
    [0 0 0 0 0
     0 0 0 0 1
     0 0 0 1 0
     0 0 0 1 1
     0 0 1 0 0
     0 0 1 0 1
     0 0 1 1 0
     0 0 1 1 1
     0 1 0 0 0
     0 1 0 0 1
     0 1 0 1 0
     0 1 0 1 1
     0 1 1 0 0
     0 1 1 0 1
     0 1 1 1 0
     0 1 1 1 1
     1 0 0 0 0
     1 0 0 0 1
     1 0 0 1 0
     1 0 0 1 1
     1 0 1 0 0
     1 0 1 0 1
     1 0 1 1 0
     1 0 1 1 1
     1 1 0 0 0
     1 1 0 0 1
     1 1 0 1 0
     1 1 0 1 1
     1 1 1 0 0
     1 1 1 0 1
     1 1 1 1 0
     1 1 1 1 1],
  :N => 1000
)
lsat[:R] = size(lsat[:response], 1)
lsat[:T] = size(lsat[:response], 2)

n = [lsat[:culm][1]; diff(lsat[:culm])]
idx = mapreduce(i -> fill(i, n[i]), vcat, 1:length(n))
lsat[:r] = lsat[:response][idx, :]


## Model Specification

model = Model(

  r = Stochastic(2,
    (beta, theta, alpha, N, T) ->
      UnivariateDistribution[(
        p = invlogit(beta * theta[i] - alpha[j]);
        Bernoulli(p)) for i in 1:N, j in 1:T
      ],
    false
  ),

  theta = Stochastic(1,
    () -> Normal(0, 1),
    false
  ),

  alpha = Stochastic(1,
    () -> Normal(0, 100),
    false
  ),

  a = Logical(1,
    alpha -> alpha .- mean(alpha)
  ),

  beta = Stochastic(
    () -> Truncated(Flat(), 0, Inf)
  )

)


## Initial Values
inits = [
  Dict(:r => lsat[:r], :alpha => zeros(lsat[:T]), :beta => 1,
       :theta => zeros(lsat[:N])),
  Dict(:r => lsat[:r], :alpha => ones(lsat[:T]), :beta => 2,
       :theta => zeros(lsat[:N]))
]


## Sampling Scheme
scheme = [AMWG(:alpha, 0.1),
          Slice(:beta, 1.0),
          Slice(:theta, 0.5)]
setsamplers!(model, scheme)


## MCMC Simulations
sim = mcmc(model, lsat, inits, 10000, burnin=2500, thin=2, chains=2)
describe(sim)
