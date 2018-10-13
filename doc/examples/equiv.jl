using Distributed
@everywhere using Mamba

## Data
equiv = Dict{Symbol, Any}(
  :group => [1, 1, 2, 2, 2, 1, 1, 1, 2, 2],
  :y =>
    [1.40 1.65
     1.64 1.57
     1.44 1.58
     1.36 1.68
     1.65 1.69
     1.08 1.31
     1.09 1.43
     1.25 1.44
     1.25 1.39
     1.30 1.52]
)
equiv[:N] = size(equiv[:y], 1)
equiv[:P] = size(equiv[:y], 2)

equiv[:T] = [equiv[:group] 3 - equiv[:group]]


## Model Specification
model = Model(

  y = Stochastic(2,
    (delta, mu, phi, pi, s2_1, T) ->
      begin
        sigma = sqrt(s2_1)
        UnivariateDistribution[
          begin
            m = mu + (-1)^(T[i, j] - 1) * phi / 2 + (-1)^(j - 1) * pi / 2 +
                delta[i, j]
            Normal(m, sigma)
          end
          for i in 1:10, j in 1:2
        ]
      end,
    false
  ),

  delta = Stochastic(2,
    s2_2 -> Normal(0, sqrt(s2_2)),
    false
  ),

  mu = Stochastic(
    () -> Normal(0, 1000)
  ),

  phi = Stochastic(
    () -> Normal(0, 1000)
  ),

  theta = Logical(
    phi -> exp(phi)
  ),

  pi = Stochastic(
    () -> Normal(0, 1000)
  ),

  s2_1 = Stochastic(
    () -> InverseGamma(0.001, 0.001)
  ),

  s2_2 = Stochastic(
    () -> InverseGamma(0.001, 0.001)
  ),

  equiv = Logical(
    theta -> Int(0.8 < theta < 1.2)
  )

)


## Initial Values
inits = [
  Dict(:y => equiv[:y], :delta => zeros(10, 2), :mu => 0, :phi => 0,
       :pi => 0, :s2_1 => 1, :s2_2 => 1),
  Dict(:y => equiv[:y], :delta => zeros(10, 2), :mu => 10, :phi => 10,
       :pi => 10, :s2_1 => 10, :s2_2 => 10)
]


## Sampling Scheme
scheme = [NUTS(:delta),
          Slice([:mu, :phi, :pi], 1.0),
          Slice([:s2_1, :s2_2], 1.0, Univariate)]
setsamplers!(model, scheme)


## MCMC Simulations
sim = mcmc(model, equiv, inits, 12500, burnin=2500, thin=2, chains=2)
describe(sim)
