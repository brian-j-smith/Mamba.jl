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


#################### Discrete Support Grids ####################

typealias GridUnivariateDistribution
          Union(Bernoulli, Binomial, Categorical, DiscreteUniform,
                Hypergeometric, NoncentralHypergeometric)

grid(d::GridUnivariateDistribution) = Float64[minimum(d):maximum(d)]
grid(d::Distribution) =
  error("discrete grid not available for ", typeof(d), " distributions")


#################### PDMatDistribution ####################

typealias PDMatDistribution Union(InverseWishart, Wishart)

function link(D::PDMatDistribution, X::Matrix, transform::Bool=true)
  n = dim(D)
  value = similar(X, Int(n * (n + 1) / 2))
  k = 1
  if transform
    U = chol(X)
    for i in 1:n
      value[k] = log(U[i,i])
      k += 1
    end
    for i in 1:n, j in (i+1):n
      value[k] = U[i,j]
      k += 1
    end
  else
    for i in 1:n, j in i:n
      value[k] = X[i,j]
      k += 1
    end
  end
  value
end

function invlink(D::PDMatDistribution, x::Vector, transform::Bool=true)
  n = dim(D)
  value = zeros(VariateType, n, n)
  k = 1
  if transform
    for i in 1:n
      value[i,i] = exp(x[k])
      k += 1
    end
    for i in 1:n, j in (i+1):n
      value[i,j] = x[k]
      k += 1
    end
    return At_mul_B(value, value)
  else
    for i in 1:n, j in i:n
      value[i,j] = value[j,i] = x[k]
      k += 1
    end
    return value
  end
end

function logpdf{T<:Real}(D::PDMatDistribution, X::Matrix{T}, transform::Bool)
  value = logpdf(D, X)
  if transform && isfinite(value)
    U = chol(X)
    n = dim(D)
    for i in 1:n
      u = U[i,i]
      value += log(2 * u) + (n - i + 1) * log(u)
    end
  end
  value
end


#################### TransformDistribution ####################

typealias TransformDistribution{T<:ContinuousUnivariateDistribution}
  Union(T, Truncated{T})

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
