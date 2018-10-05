using Mamba

@everywhere extensions = quote
  using Distributions
  import Distributions: location, scale, skewness, kurtosis, minimum, maximum,
         quantile

  struct GK <: ContinuousUnivariateDistribution
    A::Float64
    B::Float64
    g::Float64
    k::Float64
    c::Float64

    function GK(A::Real, B::Real, g::Real, k::Real, c::Real)
      ## check args
      0.0 <= c < 1.0 || throw(ArgumentError("c must be in [0, 1)"))

      ## create distribution
      new(A, B, g, k, c)
    end

    GK(A::Real, B::Real, g::Real, k::Real) = GK(A, B, g, k, 0.8)
  end

  ## Parameters
  location(d::GK) = d.A
  scale(d::GK) = d.B
  skewness(d::GK) = d.g
  kurtosis(d::GK) = d.k
  asymmetry(d::GK) = d.c

  minimum(d::GK) = -Inf
  maximum(d::GK) = Inf

  function quantile(d::GK, p::Float64)
    z = quantile(Normal(), p)
    z2gk(d, z)
  end

  function z2gk(d::GK, z::Float64)
    term1 = exp(-skewness(d) * z)
    term2 = (1.0 + asymmetry(d) * (1.0 - term1) / (1.0 + term1))
    term3 = (1.0 + z^2)^kurtosis(d)
    location(d) + scale(d) * z * term2 * term3
  end
end

@everywhere eval(extensions)

d = GK(3, 1, 2, 0.5)
x = rand(d, 1000)

allingham = Dict{Symbol, Any}(
  :x => x
)

model = Model(
  x = Stochastic(1, (A, B, g, k) -> GK(A, B, g, k), false),
  A = Stochastic(() -> Uniform(0, 10)),
  B = Stochastic(() -> Uniform(0, 10)),
  g = Stochastic(() -> Uniform(0, 10)),
  k = Stochastic(() -> Uniform(0, 10))
)

inits = [
  Dict{Symbol, Any}(:x => x, :A => 3.5, :B => 0.5,
                    :g => 2.0, :k => 0.5),
  Dict{Symbol, Any}(:x => x, :A => median(x),
                    :B => sqrt(var(x)),
                    :g => 1.0, :k => 1.0),
  Dict{Symbol, Any}(:x => x, :A => median(x),
                    :B => diff(quantile(x, [0.25, 0.75]))[1],
                    :g => mean((x - mean(x)).^3) / (var(x)^(3 / 2)),
                    :k => rand())
]

sigma1 = 0.05
sigma2 = 0.5
stats = x -> quantile(x, [0.1, 0.25, 0.5, 0.75, 0.9])
epsilon = 0.1

scheme = [ABC([:A, :B, :k],
              sigma1, stats, epsilon, maxdraw=50, decay=0.75, randeps=true),
          ABC(:g, sigma2, stats, epsilon, maxdraw=50, decay=0.75)]
setsamplers!(model, scheme)

sim = mcmc(model, allingham, inits, 10000, burnin=2500, chains=3)
describe(sim)
p = plot(sim)
draw(p, filename="gk")
