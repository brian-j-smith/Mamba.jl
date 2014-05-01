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


#################### TransformDistribution ####################

typealias TransformDistribution{T<:ContinuousUnivariateDistribution}
  Union(T, Truncated{T, Continuous})

function minimum(d::Truncated)
  max(d.lower, minimum(d.untruncated))
end

function maximum(d::Truncated)
  min(d.upper, maximum(d.untruncated))
end

function link(d::TransformDistribution, x)
  a, b = minimum(d), maximum(d)
  if a > -Inf && b < Inf
    logit((x - a) ./ (b - a))
  elseif a > -Inf
    log(x - a)
  elseif b < Inf
    log(b - x)
  else
    x
  end
end

function invlink(d::TransformDistribution, x)
  a, b = minimum(d), maximum(d)
  if a > -Inf && b < Inf
    (b - a) * invlogit(x) + a
  elseif a > -Inf
    exp(x) + a
  elseif b < Inf
    b - exp(x)
  else
    x
  end
end

function logpdf(d::TransformDistribution, x::Real, transform::Bool)
  value = logpdf(d, x)
  if transform
    a, b = minimum(d), maximum(d)
    if a > -Inf && b < Inf
      y = (x - a) / (b - x)
      value += log((b - a) * y / (y + 1.0)^2)
    elseif a > -Inf
      value += log(x - a)
    elseif b < Inf
      value += log(b - x)
    end
  end
  value
end
