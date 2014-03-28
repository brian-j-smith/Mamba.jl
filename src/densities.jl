#################### Fallback Methods ####################

link(d::Distribution, x) = x
invlink(d::Distribution, x) = x
logpdf(d::UnivariateDistribution, x, transform::Bool) = logpdf(d, x)
logpdf(d::MultivariateDistribution, x, transform::Bool) = logpdf(d, x)


#################### Null Distribution ####################

type NullDistribution <: Distribution end


#################### Flat Distribution ####################

immutable Flat <: Distribution
  length::Integer
  Flat(length::Real) = new(integer(length))
end

insupport{T<:Real}(d::Flat, x::Union(T, Vector{T})) = d.length == length(x)

function logpdf{T<:Real}(d::Flat, x::Union(T, Vector{T}), transform::Bool)
  d.length == length(x) ? 0.0 : throw(BoundsError())
end


#################### PositiveDistribution ####################

typealias PositiveDistribution
  Union(BetaPrime, Chi, Chisq, Erlang, Exponential, FDist, Gamma, InverseGamma,
        InverseGaussian, Kolmogorov, LogNormal, NoncentralChisq, NoncentralF,
        Rayleigh, Weibull)

link(d::PositiveDistribution, x) = log(x)

invlink(d::PositiveDistribution, x) = exp(x)

function logpdf(d::PositiveDistribution, x::Real, transform::Bool)
  transform ? logpdf(d, x) + log(x) : logpdf(d, x)
end


#################### RightDistribution ####################

typealias RightDistribution Union(Levy, Pareto)

link(d::RightDistribution, x) = log(x - minimum(d))

invlink(d::RightDistribution, x) = exp(x) + minimum(d)

function logpdf(d::RightDistribution, x::Real, transform::Bool)
  transform ? logpdf(d, x) + log(x - minimum(d)) : logpdf(x)
end


#################### UnitDistribution ####################

typealias UnitDistribution
  Union(Arcsine, Beta, Cosine, KSOneSided, NoncentralBeta)

link(d::UnitDistribution, x) = logit(x)

invlink(d::UnitDistribution, x) = invlogit(x)

function logpdf(d::UnitDistribution, x::Real, transform::Bool)
  if transform
    y = x / (1.0 - x)
    logpdf(d, x) + log(y / (y + 1.0)^2)
  else
    logpdf(d, x)
  end
end


#################### RangeDistribution ####################

typealias RangeDistribution Union(KSDist, TriangularDist, Uniform)

function link(d::RangeDistribution, x)
  a, b = minimum(d), maximum(d)
  logit((x - a) / (b - a))
end

function invlink(d::RangeDistribution, x)
  a, b = minimum(d), maximum(d)
  (b - a) * invlogit(x) + a
end

function logpdf(d::RangeDistribution, x::Real, transform::Bool)
  if transform
    a, b = minimum(d), maximum(d)
    y = (x - a) / (b - x)
    logpdf(d, x) + log((b - a) * y / (y + 1.0)^2)
  else
    logpdf(d, x)
  end
end
