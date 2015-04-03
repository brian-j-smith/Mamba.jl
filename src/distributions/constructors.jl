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

@distr_fallback :(Bernoulli(Float64(p)))


#################### Beta ####################

@distr_fallback :(Beta(Float64(a), Float64(b)))
@distr_fallback :(Beta(Float64(a)))


#################### BetaPrime ####################

@distr_fallback :(BetaPrime(Float64(a), Float64(b)))


#################### Binomial ####################

@distr_fallback :(Binomial(Int(n), Float64(p)))
@distr_fallback :(Binomial(Int(size)))


#################### Categorical ####################

Distributions.Categorical(p::VectorVariate) = Categorical(p.value)
Distributions.Categorical(k::UniVariate) = Categorical(Int(k))


#################### Cauchy ####################

@distr_fallback :(Cauchy(Float64(location), Float64(scale)))
@distr_fallback :(Cauchy(Float64(location)))


#################### Chi ####################

@distr_fallback :(Chi(Float64(df)))


#################### Chisq ####################

@distr_fallback :(Chisq(Float64(df)))


#################### Cosine ####################

## No constructor arguments


#################### DiscreteUniform ####################

@distr_fallback :(DiscreteUniform(Int(a), Int(b)))
@distr_fallback :(DiscreteUniform(Int(b)))


#################### Edgeworth ####################

Distributions.EdgeworthZ(d::UnivariateDistribution, n) =
  EdgeworthZ(d, Float64(n))


#################### Erlang ####################

@distr_fallback :(Erlang(Float64(shape), Float64(scale)))
@distr_fallback :(Erlang(Float64(shape)))


#################### Exponential ####################

@distr_fallback :(Exponential(Float64(scale)))


#################### FDist ####################

@distr_fallback :(FDist(Float64(d1), Float64(d2)))


#################### Gamma ####################

@distr_fallback :(Gamma(Float64(shape), Float64(scale)))
@distr_fallback :(Gamma(Float64(shape)))


#################### Geometric ####################

@distr_fallback :(Geometric(Float64(p)))


#################### Gumbel ####################

@distr_fallback :(Gumbel(Float64(mu), Float64(beta)))


#################### Hypergeometric ####################

@distr_fallback :(Hypergeometric(Float64(s), Float64(f), Float64(n)))


#################### InverseGamma ####################

@distr_fallback :(InverseGamma(Float64(shape), Float64(scale)))


#################### InverseGaussian ####################

@distr_fallback :(InverseGaussian(Float64(mu), Float64(lambda)))


#################### Kolmogorov ####################

## No constructor arguments


#################### KSDist ####################

@distr_fallback :(KSDist(Int(n)))


#################### KSOneSided ####################

@distr_fallback :(KSOneSided(Int(n)))


#################### Laplace ####################

@distr_fallback :(Laplace(Float64(location), Float64(scale)))
@distr_fallback :(Laplace(Float64(location)))


#################### Levy ####################

@distr_fallback :(Levy(Float64(location), Float64(scale)))
@distr_fallback :(Levy(Float64(location)))


#################### Logistic ####################

@distr_fallback :(Logistic(Float64(location), Float64(scale)))
@distr_fallback :(Logistic(Float64(location)))


#################### LogNormal ####################

@distr_fallback :(LogNormal(Float64(ml), Float64(sdl)))
@distr_fallback :(LogNormal(Float64(ml)))


#################### NegativeBinomial ####################

@distr_fallback :(NegativeBinomial(Float64(r), Float64(p)))


#################### NoncentralBeta ####################

@distr_fallback :(NoncentralBeta(Float64(a), Float64(b), Float64(nc)))


#################### NoncentralChisq ####################

@distr_fallback :(NoncentralChisq(Float64(df), Float64(nc)))


#################### NoncentralF ####################

@distr_fallback :(NoncentralF(Float64(n), Float64(d), Float64(nc)))


#################### NoncentralT ####################

@distr_fallback :(NoncentralT(Float64(df), Float64(nc)))


#################### Normal ####################

@distr_fallback :(Normal(Float64(μ), Float64(σ)))
@distr_fallback :(Normal(Float64(μ)))


#################### NormalCanon ####################

@distr_fallback :(NormalCanon(Float64(h), Float64(prec)))


#################### Pareto ####################

@distr_fallback :(Pareto(Float64(scale), Float64(shape)))
@distr_fallback :(Pareto(Float64(scale)))


#################### Poisson ####################

@distr_fallback :(Poisson(Float64(lambda)))


#################### Rayleigh ####################

@distr_fallback :(Rayleigh(Float64(scale)))


#################### Skellam ####################

@distr_fallback :(Skellam(Float64(mu1), Float64(mu2)))


#################### TDist ####################

@distr_fallback :(TDist(Float64(df)))


#################### TriangularDist ####################

@distr_fallback :(TriangularDist(Float64(location), Float64(scale)))
@distr_fallback :(TriangularDist(Float64(location)))


#################### Uniform ####################

@distr_fallback :(Uniform(Float64(a), Float64(b)))


#################### Weibull ####################

@distr_fallback :(Weibull(Float64(shape), Float64(scale)))
@distr_fallback :(Weibull(Float64(shape)))


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
