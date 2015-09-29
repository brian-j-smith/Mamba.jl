#################### DistributionStruct ####################

#################### Base Methods ####################

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
