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

Distributions.Categorical(p::VariateVector) = Categorical(p.value)
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

#################### GenericMvNormal ####################

Distributions.MvNormal(μ, C::PDMat) = MvNormal(convert(Vector{Float64}, μ), C)
Distributions.MvNormal(Σ) = MvNormal(convert(Matrix{Float64}, Σ))
Distributions.MvNormal(μ, Σ) =
  MvNormal(convert(Vector{Float64}, μ), convert(Matrix{Float64}, Σ))

Distributions.DiagNormal(μ, C::PDiagMat) =
  DiagNormal(convert(Vector{Float64}, μ), C)
Distributions.DiagNormal(μ, σ) =
  DiagNormal(convert(Vector{Float64}, μ), convert(Vector{Float64}, σ))

Distributions.IsoNormal(d::Integer, σ::UniVariate) = IsoNormal(d, σ.value)
Distributions.IsoNormal(μ, C::ScalMat) =
  IsoNormal(convert(Vector{Float64}, μ), C)
Distributions.IsoNormal(μ, σ) =
  IsoNormal(convert(Vector{Float64}, μ), convert(Float64, σ))


#################### Dirichlet ####################

Distributions.Dirichlet(d::Integer, alpha::UniVariate) =
  Dirichlet(d, alpha.value)
Distributions.Dirichlet(alpha) = Dirichlet(convert(Vector{Float64}, alpha))


#################### Multinomial ####################

Distributions.Multinomial(n::Integer, p) =
  Multinomial(n, convert(Vector{Float64}, p))


#################### GenericMvNormalCanon ####################

Distributions.MvNormalCanon(h, J::PDMat) =
  MvNormalCanon(convert(Vector{Float64}, h), J)
Distributions.MvNormalCanon(J) = MvNormalCanon(convert(Matrix{Float64}, J))
Distributions.MvNormalCanon(h, J) =
  MvNormalCanon(convert(Vector{Float64}, h), convert(Matrix{Float64}, J))

Distributions.DiagNormalCanon(h, J::PDiagMat) =
  DiagNormalCanon(convert(Vector{Float64}, h), J)
Distributions.DiagNormalCanon(J) = DiagNormalCanon(convert(Vector{Float64}, J))
Distributions.DiagNormalCanon(h, J) =
  DiagNormalCanon(convert(Vector{Float64}, h), convert(Vector{Float64}, J))

Distributions.IsoNormalCanon(d::Integer, prec::UniVariate) =
  IsoNormalCanon(d, prec.value)
Distributions.IsoNormalCanon(h, J::ScalMat) =
  IsoNormalCanon(convert(Vector{Float64}, h), J)
Distributions.IsoNormalCanon(h, prec) =
  IsoNormalCanon(convert(Vector{Float64}, h), convert(Float64, prec))


#################### GenericMvTDist ####################

Distributions.MvTDist(df, μ, C::PDMat) =
  MvTDist(convert(Float64, df), convert(Vector{Float64}, μ), C)
Distributions.MvTDist(df, C::PDMat) = MvTDist(convert(Float64, df), C)
Distributions.MvTDist(df, μ, Σ) =
  MvTDist(convert(Float64, df), convert(Vector{Float64}, μ),
    convert(Matrix{Float64}, Σ))
Distributions.MvTDist(df, Σ) =
  MvTDist(convert(Float64, df), convert(Matrix{Float64}, Σ))

Distributions.DiagTDist(df, μ, C::PDiagMat) =
  DiagTDist(convert(Float64, df), convert(Vector{Float64}, μ), C)
Distributions.DiagTDist(df, C::PDiagMat) =
  DiagTDist(conert(Float64, df), C::PDiagMat)
Distributions.DiagTDist(df, μ, σ) =
  DiagTDist(convert(Float64, df), convert(Vector{Float64}, μ),
            convert(Vector{Float64}, σ))

Distributions.IsoTDist(df, d::Integer, σ::Union(Real, UniVariate)) =
  IsoTDist(convert(Float64, df), d, convert(Float64, σ))
Distributions.IsoTDist(df, μ, C::ScalMat) =
  IsoTDist(convert(Float64, df), convert(Vector{Float64}, μ), C)
Distributions.IsoTDist(df, C::ScalMat) =
  IsoTDist(convert(Float64, df), C::ScalMat)
Distributions.IsoTDist(df, μ, σ) =
  IsoTDist(convert(Float64, df), convert(Vector{Float64}, μ),
           convert(Float64, σ))


#################### VonMisesFisher ####################

Distributions.VonMisesFisher(mu, kappa) =
  VonMisesFisher(convert(Vector{Float64}, mu), convert(Float64, kappa))


######################################################################
# Matrix Distributions
######################################################################

#################### InverseWishart ####################

Distributions.InverseWishart(nu, Psi) =
  InverseWishart(convert(Float64, nu), convert(Matrix{Float64}, Psi))


#################### Wishart ####################

Distributions.Wishart(nu, S) =
  Wishart(convert(Float64, nu), convert(Matrix{Float64}, S))
