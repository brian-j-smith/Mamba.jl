#################### TransformDistribution ####################

typealias TransformDistribution{T<:ContinuousUnivariateDistribution}
  Union{T, Truncated{T}}

function link(d::TransformDistribution, x::Real, transform::Bool)
  y = x
  if transform
    a, b = minimum(d), maximum(d)
    lowerbounded, upperbounded = isfinite(a), isfinite(b)
    if lowerbounded && upperbounded
      y = logit((x - a) / (b - a))
    elseif lowerbounded
      y = log(x - a)
    elseif upperbounded
      y = log(b - x)
    end
  end
  y
end

function invlink(d::TransformDistribution, x::Real, transform::Bool)
  y = x
  if transform
    a, b = minimum(d), maximum(d)
    lowerbounded, upperbounded = isfinite(a), isfinite(b)
    if lowerbounded && upperbounded
      y = (b - a) * invlogit(x) + a
    elseif lowerbounded
      y = exp(x) + a
    elseif upperbounded
      y = b - exp(x)
    end
  end
  y
end

function logpdf(d::TransformDistribution, x::Real, transform::Bool)
  lp = logpdf(d, x)
  if transform
    a, b = minimum(d), maximum(d)
    lowerbounded, upperbounded = isfinite(a), isfinite(b)
    if lowerbounded && upperbounded
      y = (x - a) / (b - x)
      lp += log((b - a) * y / (y + 1.0)^2)
    elseif lowerbounded
      lp += log(x - a)
    elseif upperbounded
      lp += log(b - x)
    end
  end
  lp
end
