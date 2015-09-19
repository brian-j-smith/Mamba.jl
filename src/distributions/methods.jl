#################### DistributionStruct Methods ####################

dims(d::DistributionStruct) = size(d)

function dims(D::Array{MultivariateDistribution})
  if length(D) > 0
    n = length(D[1])
    all(i -> length(D[i]) == n, 2:length(D)) ||
      error("lengths of distribution array elements differ")
  else
    n = 0
  end
  size(D)..., n
end


#################### Link Fallbacks ####################

link(d::Distribution, x, transform::Bool=true) = x

function link(D::Array{UnivariateDistribution}, X::Array, transform::Bool=true)
  Y = similar(X)
  map!(i -> link(D[i], X[i], transform), Y, 1:length(D))
end

link(D::Array{MultivariateDistribution}, X::Array, transform::Bool=true) = X


invlink(d::Distribution, x, transform::Bool=true) = x

function invlink(D::Array{UnivariateDistribution}, X::Array,
                 transform::Bool=true)
  Y = similar(X)
  map!(i -> invlink(D[i], X[i], transform), Y, 1:length(D))
end

invlink(D::Array{MultivariateDistribution}, X::Array, transform::Bool=true) = X


#################### Logpdf Fallbacks ####################

function logpdf(d::Distribution, x, transform::Bool)
  insupport(d, x) ? logpdf(d, x) : -Inf
end

function logpdf(d::UnivariateDistribution, X::Array{Float64}, transform::Bool)
  map(x -> logpdf(d, x, transform), X)
end

function logpdf(D::Array{UnivariateDistribution}, X::Array{Float64},
                transform::Bool=false)
  Y = similar(D, Float64)
  map!(i -> logpdf(D[i], X[i], transform), Y, 1:length(D))
end

function logpdf(D::Array{MultivariateDistribution}, X::Array{Float64},
                transform::Bool=false)
  Y = similar(D, Float64)
  map!(i -> logpdf(D[i], vec(X[ind2sub(D, i)..., :]), transform), Y,
       1:length(D))
end


#################### Discrete Support Grids ####################

typealias GridUnivariateDistribution
          Union(Bernoulli, Binomial, Categorical, DiscreteUniform,
                Hypergeometric, NoncentralHypergeometric)

grid(d::GridUnivariateDistribution) =
  collect(UnitRange{Float64}(minimum(d), maximum(d)))

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
  value = zeros(Float64, n, n)
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

function logpdf(D::PDMatDistribution, X::Matrix{Float64}, transform::Bool)
  value = logpdf(D, X)
  if transform && isfinite(value)
    U = chol(X)
    n = dim(D)
    for i in 1:n
      value += (n - i + 2) * log(U[i,i])
    end
    value += n * log(2)
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

function logpdf(d::TransformDistribution, x::Float64, transform::Bool)
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
