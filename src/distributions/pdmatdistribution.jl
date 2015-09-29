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
