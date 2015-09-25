#################### DistributionStruct Methods ####################

dims(d::DistributionStruct) = size(d)

function dims(D::Array{MultivariateDistribution})
  size(D)..., mapreduce(length, max, D)
end


#################### List Fallbacks ####################

unlist_sub(d::Distribution, x) = x

unlist_sub(D::Array{UnivariateDistribution}, x) = x

function unlist_sub(D::Array{MultivariateDistribution}, X::Array)
  y = similar(X, length(X))
  offset = 0
  for sub in CartesianRange(size(D))
    n = length(D[sub])
    inds = 1:n
    y[inds + offset] = X[sub, inds]
    offset += n
  end
  resize!(y, offset)
end


relist_sub(d::Distribution, x) = x

relist_sub(D::Array{UnivariateDistribution}, x) = x

function relist_sub(D::Array{MultivariateDistribution}, x::Array)
  Y = similar(x, dims(D))
  offset = 0
  for sub in CartesianRange(size(D))
    n = length(D[sub])
    inds = 1:n
    Y[sub, inds] = x[inds + offset]
    offset += n
  end
  Y
end


#################### Link Fallbacks ####################

link_sub(d::Distribution, x, transform::Bool) = x

function link_sub(D::Array{UnivariateDistribution}, X::Array, transform::Bool)
  Y = similar(X)
  map!(i -> link_sub(D[i], X[i], transform), Y, 1:length(D))
end

link_sub(D::Array{MultivariateDistribution}, X::Array, transform::Bool) = X


invlink_sub(d::Distribution, x, transform::Bool) = x

function invlink_sub(D::Array{UnivariateDistribution}, X::Array,
                     transform::Bool)
  Y = similar(X)
  map!(i -> invlink_sub(D[i], X[i], transform), Y, 1:length(D))
end

invlink_sub(D::Array{MultivariateDistribution}, X::Array, transform::Bool) = X


#################### Logpdf Fallbacks ####################

function logpdf_sub(d::Distribution, x, transform::Bool)
  insupport(d, x) ? logpdf(d, x) : -Inf
end

function logpdf_sub(d::UnivariateDistribution, X::Array{Float64},
                    transform::Bool)
  map(x -> logpdf_sub(d, x, transform), X)
end

function logpdf_sub(D::Array{UnivariateDistribution}, X::Array{Float64},
                    transform::Bool)
  Y = similar(D, Float64)
  map!(i -> logpdf_sub(D[i], X[i], transform), Y, 1:length(D))
end

function logpdf_sub(D::Array{MultivariateDistribution}, X::Array{Float64},
                    transform::Bool)
  Y = similar(D, Float64)
  for sub in CartesianRange(size(D))
    d = D[sub]
    Y[sub] = logpdf_sub(d, vec(X[sub, 1:length(d)]), transform)
  end
  Y
end


#################### Rand Fallbacks ####################

rand_sub(d::Distribution) = rand(d)

rand_sub(D::Array{UnivariateDistribution}) = map(rand, D)

function rand_sub(D::Array{MultivariateDistribution})
  x = fill(NaN, dims(D))
  for sub in CartesianRange(size(D))
    d = D[sub]
    x[sub, 1:length(d)] = rand(d)
  end
  x
end


#################### Discrete Support Grids ####################

typealias GridUnivariateDistribution
          Union{Bernoulli, Binomial, Categorical, DiscreteUniform,
                Hypergeometric, NoncentralHypergeometric}


#################### PDMatDistribution ####################

typealias PDMatDistribution Union{InverseWishart, Wishart}

function unlist_sub(d::PDMatDistribution, X::Matrix)
  n = dim(d)
  y = similar(X, Int(n * (n + 1) / 2))
  k = 0
  for i in 1:n, j in i:n
    k += 1
    y[k] = X[i,j]
  end
  y
end

function relist_sub(d::PDMatDistribution, x::Array)
  n = dim(d)
  Y = similar(x, n, n)
  k = 0
  for i in 1:n, j in i:n
    k += 1
    Y[i,j] = Y[j,i] = x[k]
  end
  Y
end

function link_sub(d::PDMatDistribution, X::Matrix{Float64}, transform::Bool)
  if transform
    n = dim(d)
    Y = zeros(n, n)
    U = chol(X)
    for i in 1:n
      Y[i,i] = log(U[i,i])
    end
    for i in 1:n, j in (i+1):n
      Y[i,j] = U[i,j]
    end
  else
    Y = X
  end
  Y
end

function invlink_sub(d::PDMatDistribution, X::Matrix{Float64}, transform::Bool)
  if transform
    n = dim(d)
    U = zeros(n, n)
    for i in 1:n
      U[i,i] = exp(X[i,i])
    end
    for i in 1:n, j in (i+1):n
      U[i,j] = X[i,j]
    end
    Y = At_mul_B(U, U)
  else
    Y = X
  end
  Y
end

function logpdf_sub(d::PDMatDistribution, X::Matrix{Float64}, transform::Bool)
  value = logpdf(d, X)
  if transform && isfinite(value)
    U = chol(X)
    n = dim(d)
    for i in 1:n
      value += (n - i + 2) * log(U[i,i])
    end
    value += n * log(2)
  end
  value
end


#################### TransformDistribution ####################

typealias TransformDistribution{T<:ContinuousUnivariateDistribution}
  Union{T, Truncated{T}}

function link_sub(d::TransformDistribution, x, transform::Bool)
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

function invlink_sub(d::TransformDistribution, x, transform::Bool)
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

function logpdf_sub(d::TransformDistribution, x::Float64, transform::Bool)
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
