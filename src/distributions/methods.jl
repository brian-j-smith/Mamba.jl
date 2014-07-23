#################### Fallback Methods ####################

link(d::Distribution, x, transform::Bool=true) = x
invlink(d::Distribution, x, transform::Bool=true) = x

function logpdf(d::UnivariateDistribution, x, transform::Bool)
  all(insupport(d, x)) ? logpdf(d, x) : -Inf
end

function logpdf(d::MultivariateDistribution, x, transform::Bool)
  all(insupport(d, x)) ? logpdf(d, x) : -Inf
end

function logpdf(d::MatrixDistribution, x, transform::Bool)
  all(insupport(d, x)) ? logpdf(d, x) : -Inf
end


#################### Truncated ####################

function rand(d::Truncated)
  if d.nc > 0.25
    while true
      r = rand(d.untruncated)
      if d.lower <= r <= d.upper
        return r
      end
    end
  else
    return quantile(d.untruncated, cdf(d.untruncated, d.lower) + rand() * d.nc)
  end
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

function link(d::TransformDistribution, x, transform::Bool=true)
  if transform
    a, b = minimum(d), maximum(d)
    lowerbounded, upperbounded = isfinite(a), isfinite(b)
    if lowerbounded && upperbounded
      return logit((x - a) / (b - a))
    elseif lowerbounded
      return log(x - a)
    elseif upperbounded
      return log(b - x)
    else
      return x
    end
  else
     return x
  end
end

function invlink(d::TransformDistribution, x, transform::Bool=true)
  if transform
    a, b = minimum(d), maximum(d)
    lowerbounded, upperbounded = isfinite(a), isfinite(b)
    if lowerbounded && upperbounded
      return (b - a) * invlogit(x) + a
    elseif lowerbounded
      return exp(x) + a
    elseif upperbounded
      return b - exp(x)
    else
      return x
    end
  else
    return x
  end
end

function logpdf(d::TransformDistribution, x::Real, transform::Bool)
  insupport(d, x) || return -Inf
  value = logpdf(d, x)
  if transform
    a, b = minimum(d), maximum(d)
    lowerbounded, upperbounded = isfinite(a), isfinite(b)
    if lowerbounded && upperbounded
      y = (x - a) / (b - x)
      value += log((b - a) * y / (y + 1.0)^2)
    elseif lowerbounded
      value += log(x - a)
    elseif upperbounded
      value += log(b - x)
    end
  end
  value
end
