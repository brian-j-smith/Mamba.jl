#################### DistributionStruct Methods ####################

dims(d::DistributionStruct) = size(d)

function dims(D::Array{MultivariateDistribution})
  size(D)..., mapreduce(length, max, D)
end


#################### List Fallbacks ####################

unlist(d::Distribution, x) = x

unlist_sub(d::Distribution, x) = unlist(d, x)

unlist_sub(D::Array{UnivariateDistribution}, x) = x

function unlist_sub(D::Array{MultivariateDistribution}, X::DenseArray)
  y = similar(X, length(X))
  offset = 0
  for sub in CartesianRange(size(D))
    n = length(D[sub])
    inds = 1:n
    y[offset + inds] = X[sub, inds]
    offset += n
  end
  resize!(y, offset)
end


relist(d::Distribution, x) = x

relist_sub(d::Distribution, x) = relist(d, x)

relist_sub(D::Array{UnivariateDistribution}, x) = x

function relist_sub(D::Array{MultivariateDistribution}, x::DenseArray)
  Y = similar(x, dims(D))
  offset = 0
  for sub in CartesianRange(size(D))
    n = length(D[sub])
    inds = 1:n
    Y[sub, inds] = x[offset + inds]
    offset += n
  end
  Y
end


#################### Link Fallbacks ####################

link(d::Distribution, x, transform::Bool) = x

link_sub(d::Distribution, x, transform::Bool) = link(d, x, transform)

function link_sub(d::UnivariateDistribution, X::DenseArray, transform::Bool)
  Y = similar(X, Float64)
  map!(x -> link_sub(d, x, transform), Y, X)
end

function link_sub(D::Array{UnivariateDistribution}, X::DenseArray,
                  transform::Bool)
  Y = similar(X, Float64)
  map!(i -> link_sub(D[i], X[i], transform), Y, 1:length(D))
end

function link_sub(D::Array{MultivariateDistribution}, X::DenseArray,
                  transform::Bool)
  Y = similar(X, Float64)
  for sub in CartesianRange(size(D))
    d = D[sub]
    inds = 1:length(d)
    Y[sub, inds] = link_sub(d, X[sub, inds], transform)
  end
  Y
end


invlink(d::Distribution, x, transform::Bool) = x

invlink_sub(d::Distribution, x, transform::Bool) = invlink(d, x, transform)

function invlink_sub(d::UnivariateDistribution, X::DenseArray, transform::Bool)
  Y = similar(X, Float64)
  map!(x -> invlink_sub(d, x, transform), Y, X)
end

function invlink_sub(D::Array{UnivariateDistribution}, X::DenseArray,
                     transform::Bool)
  Y = similar(X, Float64)
  map!(i -> invlink_sub(D[i], X[i], transform), Y, 1:length(D))
end

function invlink_sub(D::Array{MultivariateDistribution}, X::DenseArray,
                     transform::Bool)
  Y = similar(X, Float64)
  for sub in CartesianRange(size(D))
    d = D[sub]
    inds = 1:length(d)
    Y[sub, inds] = invlink_sub(d, X[sub, inds], transform)
  end
  Y
end


#################### Logpdf Fallbacks ####################

logpdf(d::Distribution, x, transform::Bool) = logpdf(d, x)

function logpdf_sub(d::Distribution, x, transform::Bool)
  insupport(d, x) ? logpdf(d, x, transform) : -Inf
end

function logpdf_sub(d::UnivariateDistribution, X::DenseArray, transform::Bool)
  lp = 0.0
  for x in X
    lp += logpdf_sub(d, x, transform)
  end
  lp
end

function logpdf_sub(D::Array{UnivariateDistribution}, X::DenseArray,
                    transform::Bool)
  lp = 0.0
  for i in 1:length(D)
    lp += logpdf_sub(D[i], X[i], transform)
  end
  lp
end

function logpdf_sub(D::Array{MultivariateDistribution}, X::DenseArray,
                    transform::Bool)
  lp = 0.0
  for sub in CartesianRange(size(D))
    d = D[sub]
    lp += logpdf_sub(d, vec(X[sub, 1:length(d)]), transform)
  end
  lp
end


#################### Rand Fallbacks ####################

rand_sub(d::Distribution) = rand(d)

function rand_sub(D::Array{UnivariateDistribution})
  X = similar(D, Float64)
  map!(rand, X, D)
end

function rand_sub(D::Array{MultivariateDistribution})
  X = fill(NaN, dims(D))
  for sub in CartesianRange(size(D))
    d = D[sub]
    X[sub, 1:length(d)] = rand(d)
  end
  X
end


#################### Discrete Support Grids ####################

typealias GridUnivariateDistribution
          Union{Bernoulli, Binomial, Categorical, DiscreteUniform,
                Hypergeometric, NoncentralHypergeometric}


#################### PDMatDistribution ####################

typealias PDMatDistribution Union{InverseWishart, Wishart}

function unlist(d::PDMatDistribution, X::DenseMatrix)
  n = dim(d)
  y = similar(X, Int(n * (n + 1) / 2))
  k = 0
  for i in 1:n, j in i:n
    k += 1
    y[k] = X[i,j]
  end
  y
end

function relist(d::PDMatDistribution, x::DenseArray)
  n = dim(d)
  Y = similar(x, n, n)
  k = 0
  for i in 1:n, j in i:n
    k += 1
    Y[i,j] = Y[j,i] = x[k]
  end
  Y
end

function link(d::PDMatDistribution, X::DenseMatrix, transform::Bool)
  Y = X
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
  end
  Y
end

function invlink(d::PDMatDistribution, X::DenseMatrix, transform::Bool)
  Y = X
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
  end
  Y
end

function logpdf(d::PDMatDistribution, X::DenseMatrix, transform::Bool)
  lp = logpdf(d, X)
  if transform && isfinite(lp)
    U = chol(X)
    n = dim(d)
    for i in 1:n
      lp += (n - i + 2) * log(U[i,i])
    end
    lp += n * log(2)
  end
  lp
end


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
