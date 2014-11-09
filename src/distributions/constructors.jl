######################################################################
# Uniform Distributions
######################################################################

#################### Fallback Macro ####################

macro distr_fallback(g)
  quote
    local f = deepcopy($g)
    for i in 2:length(($g).args)
      f.args[i] = ($g).args[i].args[end]
    end
    eval(Distributions, Expr(:(=), f, $g))
  end
end


#################### Arcsine ####################

## No constructor arguments


#################### Bernoulli ####################

@distr_fallback :(Bernoulli(float64(p)))


#################### Beta ####################

@distr_fallback :(Beta(float64(a), float64(b)))
@distr_fallback :(Beta(float64(a)))


#################### BetaPrime ####################

@distr_fallback :(BetaPrime(float64(a), float64(b)))


#################### Binomial ####################

@distr_fallback :(Binomial(int(n), float64(p)))
@distr_fallback :(Binomial(int(size)))


#################### Categorical ####################

Distributions.Categorical(p::VectorVariate) = Categorical(p.value)
Distributions.Categorical(k::UniVariate) = Categorical(int(k))


#################### Cauchy ####################

@distr_fallback :(Cauchy(float64(location), float64(scale)))
@distr_fallback :(Cauchy(float64(location)))


#################### Chi ####################

@distr_fallback :(Chi(float64(df)))


#################### Chisq ####################

@distr_fallback :(Chisq(float64(df)))


#################### Cosine ####################

## No constructor arguments


#################### DiscreteUniform ####################

@distr_fallback :(DiscreteUniform(int(a), int(b)))
@distr_fallback :(DiscreteUniform(int(b)))


#################### Edgeworth ####################

Distributions.EdgeworthZ(d::UnivariateDistribution, n) =
  EdgeworthZ(d, float64(n))


#################### Erlang ####################

@distr_fallback :(Erlang(float64(shape), float64(scale)))
@distr_fallback :(Erlang(float64(shape)))


#################### Exponential ####################

@distr_fallback :(Exponential(float64(scale)))


#################### FDist ####################

@distr_fallback :(FDist(float64(d1), float64(d2)))


#################### Gamma ####################

@distr_fallback :(Gamma(float64(shape), float64(scale)))
@distr_fallback :(Gamma(float64(shape)))


#################### Geometric ####################

@distr_fallback :(Geometric(float64(p)))


#################### Gumbel ####################

@distr_fallback :(Gumbel(float64(mu), float64(beta)))


#################### Hypergeometric ####################

@distr_fallback :(Hypergeometric(float64(s), float64(f), float64(n)))


#################### InverseGamma ####################

@distr_fallback :(InverseGamma(float64(shape), float64(scale)))


#################### InverseGaussian ####################

@distr_fallback :(InverseGaussian(float64(mu), float64(lambda)))


#################### Kolmogorov ####################

## No constructor arguments


#################### KSDist ####################

@distr_fallback :(KSDist(int(n)))


#################### KSOneSided ####################

@distr_fallback :(KSOneSided(int(n)))


#################### Laplace ####################

@distr_fallback :(Laplace(float64(location), float64(scale)))
@distr_fallback :(Laplace(float64(location)))


#################### Levy ####################

@distr_fallback :(Levy(float64(location), float64(scale)))
@distr_fallback :(Levy(float64(location)))


#################### Logistic ####################

@distr_fallback :(Logistic(float64(location), float64(scale)))
@distr_fallback :(Logistic(float64(location)))


#################### LogNormal ####################

@distr_fallback :(LogNormal(float64(ml), float64(sdl)))
@distr_fallback :(LogNormal(float64(ml)))


#################### NegativeBinomial ####################

@distr_fallback :(NegativeBinomial(float64(r), float64(p)))


#################### NoncentralBeta ####################

@distr_fallback :(NoncentralBeta(float64(a), float64(b), float64(nc)))


#################### NoncentralChisq ####################

@distr_fallback :(NoncentralChisq(float64(df), float64(nc)))


#################### NoncentralF ####################

@distr_fallback :(NoncentralF(float64(n), float64(d), float64(nc)))


#################### NoncentralT ####################

@distr_fallback :(NoncentralT(float64(df), float64(nc)))


#################### Normal ####################

@distr_fallback :(Normal(float64(μ), float64(σ)))
@distr_fallback :(Normal(float64(μ)))


#################### NormalCanon ####################

@distr_fallback :(NormalCanon(float64(h), float64(prec)))


#################### Pareto ####################

@distr_fallback :(Pareto(float64(scale), float64(shape)))
@distr_fallback :(Pareto(float64(scale)))


#################### Poisson ####################

@distr_fallback :(Poisson(float64(lambda)))


#################### Rayleigh ####################

@distr_fallback :(Rayleigh(float64(scale)))


#################### Skellam ####################

@distr_fallback :(Skellam(float64(mu1), float64(mu2)))


#################### TDist ####################

@distr_fallback :(TDist(float64(df)))


#################### TriangularDist ####################

@distr_fallback :(TriangularDist(float64(location), float64(scale)))
@distr_fallback :(TriangularDist(float64(location)))


#################### Uniform ####################

@distr_fallback :(Uniform(float64(a), float64(b)))


#################### Weibull ####################

@distr_fallback :(Weibull(float64(shape), float64(scale)))
@distr_fallback :(Weibull(float64(shape)))


######################################################################
# Truncated Distributions
######################################################################

Distributions.Truncated{S<:ValueSupport}(d::UnivariateDistribution{S}, l, u) =
  Truncated(d, convert(Float64, l), convert(Float64, u))


######################################################################
# Multivariate Distributions
######################################################################

#################### MvNormal ####################

Distributions.MvNormal(μ::VectorVariate, Σ::Matrix{Float64}) =
  MvNormal(convert(Vector{Float64}, μ), Σ)
Distributions.MvNormal(μ::Vector{Float64}, Σ::MatrixVariate) =
  MvNormal(μ, convert(Matrix{Float64}, Σ))
Distributions.MvNormal(μ::VectorVariate, Σ::MatrixVariate) =
  MvNormal(convert(Vector{Float64}, μ), convert(Matrix{Float64}, Σ))

Distributions.MvNormal(μ::VectorVariate, σ::Vector{Float64}) =
  MvNormal(convert(Vector{Float64}, μ), σ)
Distributions.MvNormal(μ::Vector{Float64}, σ::VectorVariate) =
  MvNormal(μ, convert(Vector{Float64}, σ))
Distributions.MvNormal(μ::VectorVariate, σ::VectorVariate) =
  MvNormal(convert(Vector{Float64}, μ), convert(Vector{Float64}, σ))

Distributions.MvNormal(μ::VectorVariate, σ::Real) =
  MvNormal(convert(Vector{Float64}, μ), σ)
Distributions.MvNormal(μ::Vector{Float64}, σ::UniVariate) =
  MvNormal(μ, convert(Float64, σ))
Distributions.MvNormal(μ::VectorVariate, σ::UniVariate) =
  MvNormal(convert(Vector{Float64}, μ), convert(Float64, σ))

Distributions.MvNormal(Σ::MatrixVariate) =
  MvNormal(convert(Matrix{Float64}, Σ))

Distributions.MvNormal(σ::VectorVariate) =
  MvNormal(convert(Vector{Float64}, σ))

Distributions.MvNormal(d::Int, σ::UniVariate) =
  MvNormal(d, convert(Float64, σ))


#################### BDiagNormal ####################

BDiagNormal(μ::VectorVariate, Σ::Matrix{Float64}) =
  BDiagNormal(convert(Vector{Float64}, μ), Σ)
BDiagNormal(μ::Vector{Float64}, Σ::MatrixVariate) =
  BDiagNormal(μ, convert(Matrix{Float64}, Σ))
BDiagNormal(μ::VectorVariate, Σ::MatrixVariate) =
  BDiagNormal(convert(Vector{Float64}, μ), convert(Matrix{Float64}, Σ))

BDiagNormal(μ, Σ::Vector) =
  BDiagNormal(convert(Vector{Float64}, μ), Matrix{Float64}[Σ...])


#################### Dirichlet ####################

Distributions.Dirichlet(d::Integer, alpha::UniVariate) =
  Dirichlet(d, alpha.value)
Distributions.Dirichlet(alpha) = Dirichlet(convert(Vector{Float64}, alpha))


#################### Multinomial ####################

Distributions.Multinomial(n::Integer, p) =
  Multinomial(n, convert(Vector{Float64}, p))


#################### MvNormalCanon ####################

Distributions.MvNormalCanon(h::VectorVariate, J::Matrix{Float64}) =
  MvNormalCanon(convert(Vector{Float64}, h), J)
Distributions.MvNormalCanon(h::Vector{Float64}, J::MatrixVariate) =
  MvNormalCanon(h, convert(Matrix{Float64}, J))
Distributions.MvNormalCanon(h::VectorVariate, J::MatrixVariate) =
  MvNormalCanon(convert(Vector{Float64}, h), convert(Matrix{Float64}, J))

Distributions.MvNormalCanon(h::VectorVariate, prec::Vector{Float64}) =
  MvNormalCanon(convert(Vector{Float64}, h), prec)
Distributions.MvNormalCanon(h::Vector{Float64}, prec::VectorVariate) =
  MvNormalCanon(h, convert(Vector{Float64}, prec))
Distributions.MvNormalCanon(h::VectorVariate, prec::VectorVariate) =
  MvNormalCanon(convert(Vector{Float64}, h), convert(Vector{Float64}, prec))

Distributions.MvNormalCanon(h::VectorVariate, prec::Float64) =
  MvNormalCanon(convert(Vector{Float64}, h), prec)
Distributions.MvNormalCanon(h::Vector{Float64}, prec::UniVariate) =
  MvNormalCanon(h, convert(Float64, prec))
Distributions.MvNormalCanon(h::VectorVariate, prec::UniVariate) =
  MvNormalCanon(convert(Vector{Float64}, h), convert(Float64, prec))

Distributions.MvNormalCanon(J::MatrixVariate) =
  MvNormalCanon(convert(Matrix{Float64}, J))

Distributions.MvNormalCanon(prec::VectorVariate) =
  MvNormalCanon(convert(Vector{Float64}, prec))

Distributions.MvNormalCanon(d::Int, prec::UniVariate) =
  MvNormalCanon(d, convert(Float64, prec))


#################### MvTDist ####################

Distributions.MvTDist(df::UniVariate, μ::Vector{Float64}, C::PDMat) =
  MvTDist(convert(Float64, df), μ, C)
Distributions.MvTDist(df::Float64, μ::VectorVariate, C::PDMat) =
  MvTDist(df, convert(Vector{Float64}, μ), C)
Distributions.MvTDist(df::UniVariate, μ::VectorVariate, C::PDMat) =
  MvTDist(convert(Float64, df), convert(Vector{Float64}, C))

Distributions.MvTDist(df::UniVariate, C::PDMat) =
  MvTDist(convert(Float64, df), C)

Distributions.MvTDist(df::UniVariate, μ::Vector{Float64}, Σ::Matrix{Float64}) =
  MvTDist(convert(Float64, df), μ, Σ)
Distributions.MvTDist(df::Float64, μ::VectorVariate, Σ::Matrix{Float64}) =
  MvTDist(df, convert(Vector{Float64}, μ), Σ)
Distributions.MvTDist(df::Float64, μ::Vector{Float64}, Σ::MatrixVariate) =
  MvTDist(df, μ, convert(Matrix{Float64}, Σ))
Distributions.MvTDist(df::UniVariate, μ::VectorVariate, Σ::Matrix{Float64}) =
  MvTDist(convert(Float64, df), convert(Vector{Float64}, μ), Σ)
Distributions.MvTDist(df::UniVariate, μ::Vector{Float64}, Σ::MatrixVariate) =
  MvTDist(convert(Float64, df), μ, convert(Matrix{Float64}, Σ))
Distributions.MvTDist(df::Float64, μ::VectorVariate, Σ::MatrixVariate) =
  MvTDist(df, convert(Vector{Float64}, μ), convert(Matrix{Float64}, Σ))
Distributions.MvTDist(df::UniVariate, μ::VectorVariate, Σ::MatrixVariate) =
  MvTDist(convert(Float64, df), convert(Vector{Float64}, μ), convert(Matrix{Float64}, Σ))

Distributions.MvTDist(df::UniVariate, Σ::Matrix{Float64}) =
  MvTDist(convert(Float64, df), Σ)
Distributions.MvTDist(df::Float64, Σ::MatrixVariate) =
  MvTDist(df, convert(Matrix{Float64}, Σ))
Distributions.MvTDist(df::UniVariate, Σ::MatrixVariate) =
  MvTDist(convert(Float64, df), convert(Matrix{Float64}, Σ))


#################### VonMisesFisher ####################

Distributions.VonMisesFisher(mu, kappa) =
  VonMisesFisher(convert(Vector{Float64}, mu), convert(Float64, kappa))


######################################################################
# Matrix Distributions
######################################################################

#################### InverseWishart ####################

Distributions.InverseWishart(nu::UniVariate, Ψ::Matrix{Float64}) =
  InverseWishart(convert(Float64, nu), Ψ)
Distributions.InverseWishart(nu::Real, Ψ::MatrixVariate) =
  InverseWishart(nu, convert(Matrix{Float64}, Ψ))
Distributions.InverseWishart(nu::UniVariate, Ψ::MatrixVariate) =
  InverseWishart(convert(Float64, nu), convert(Matrix{Float64}, Ψ))


#################### Wishart ####################

Distributions.Wishart(nu::UniVariate, S::Matrix{Float64}) =
  Wishart(convert(Float64, nu), S)
Distributions.Wishart(nu::Real, S::MatrixVariate) =
  Wishart(nu, convert(Matrix{Float64}, S))
Distributions.Wishart(nu::UniVariate, S::MatrixVariate) =
  Wishart(convert(Float64, nu), convert(Matrix{Float64}, S))
