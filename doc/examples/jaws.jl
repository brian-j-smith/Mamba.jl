using Distributed
@everywhere using Mamba, LinearAlgebra

## Data
jaws = Dict{Symbol, Any}(
  :Y =>
   [47.8 48.8 49.0 49.7
    46.4 47.3 47.7 48.4
    46.3 46.8 47.8 48.5
    45.1 45.3 46.1 47.2
    47.6 48.5 48.9 49.3
    52.5 53.2 53.3 53.7
    51.2 53.0 54.3 54.5
    49.8 50.0 50.3 52.7
    48.1 50.8 52.3 54.4
    45.0 47.0 47.3 48.3
    51.2 51.4 51.6 51.9
    48.5 49.2 53.0 55.5
    52.1 52.8 53.7 55.0
    48.2 48.9 49.3 49.8
    49.6 50.4 51.2 51.8
    50.7 51.7 52.7 53.3
    47.2 47.7 48.4 49.5
    53.3 54.6 55.1 55.3
    46.2 47.5 48.1 48.4
    46.3 47.6 51.3 51.8],
  :age => [8.0, 8.5, 9.0, 9.5]
)
M = jaws[:M] = size(jaws[:Y], 2)
N = jaws[:N] = size(jaws[:Y], 1)
jaws[:y] = vec(jaws[:Y])
jaws[:x] = kron(ones(jaws[:N]), jaws[:age])


## Model Specification
model = Model(

  y = Stochastic(1,
    (beta0, beta1, x, Sigma) -> BDiagNormal(beta0 .+ beta1 * x, Sigma),
    false
  ),

  beta0 = Stochastic(
    () -> Normal(0, sqrt(1000))
  ),

  beta1 = Stochastic(
    () -> Normal(0, sqrt(1000))
  ),

  Sigma = Stochastic(2,
    M -> InverseWishart(4.0, Matrix{Float64}(I, M, M))
  )

)


## Initial Values
inits = [
  Dict(:y => jaws[:y], :beta0 => 40, :beta1 => 1, :Sigma => Matrix{Float64}(I, M, M)),
  Dict(:y => jaws[:y], :beta0 => 10, :beta1 => 10, :Sigma => Matrix{Float64}(I, M, M))
]


## Sampling Scheme
scheme = [Slice([:beta0, :beta1], [10, 1]),
          AMWG(:Sigma, 0.1)]
setsamplers!(model, scheme)


## MCMC Simulations
sim = mcmc(model, jaws, inits, 10000, burnin=2500, thin=2, chains=2)
describe(sim)
