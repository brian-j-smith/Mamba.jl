#################### Posterior Statistics ####################

function autocor(c::MCMCChains; lags::Vector=[1,5,10,50], relative::Bool=true)
  if relative
    lags *= c.thin
  elseif any(lags .% c.thin .!= 0)
    error("lags do not correspond to thinning interval")
  end
  labels = map(x -> "Lag " * string(x), lags)
  vals = mapslices(x -> autocor(x, lags)', c.data, [1,2])
  ChainSummary(vals, c.names, labels, header(c))
end

function batchSE(x::Vector; size::Integer=100)
  m = div(length(x), size)
  m >= 2 || error("2 or more batches needed to compute SE")
  mbar = [mean(x[i*size+(1:size)]) for i in 0:m-1]
  sem(mbar)
end

function cor(c::MCMCChains)
  ChainSummary(cor(combine(c)), c.names, c.names, header(c))
end

function describe(c::MCMCChains; batchsize::Integer=100,
                  q::Vector=[0.025, 0.25, 0.5, 0.75, 0.975])
  println(header(c))
  print("Empirical Posterior Estimates:\n")
  show(summarystats(c, batchsize=batchsize))
  print("\nQuantiles:\n")
  show(quantile(c, q=q))
end

function dic(c::MCMCChains)
  m = c.model
  nkeys = keys(m, :terminal)
  idx = indexin(labels(m, keys(m, :block)), c)

  x0 = unlist(m)

  xbar = map(i -> mean(c.data[:,i,:]), idx)
  relist!(m, xbar)
  Dhat = -2.0 * mapreduce(key -> logpdf(m[key]), +, nkeys)
  D = -2.0 * logpdf(c, nkeys)
  p = [mean(D) - Dhat, 0.5 * var(D)]

  relist!(m, x0)

  ChainSummary(hcat(Dhat + 2 * p, p), ["pD", "pV"],
               ["DIC", "Effective Parameters"], header(c))
end

function hpd(x::Vector; alpha::Real=0.05)
  n = length(x)
  m = max(1, ceil(alpha * n))

  y = sort(x)
  a = y[1:m]
  b = y[(n - m + 1):n]
  i = sortperm(b - a)[1]

  [a[i], b[i]]
end

function hpd(c::MCMCChains; alpha::Real=0.05)
  X = combine(c)
  labels = [string(100 * alpha / 2) * "%", string(100 * (1 - alpha / 2)) * "%"]
  vals = mapreduce(i -> hpd(X[:,i], alpha=alpha), hcat, 1:size(X, 2))
  ChainSummary(vals', c.names, labels, header(c))
end

function logpdf{T<:String}(c::MCMCChains, nkeys::Vector{T})
  m = c.model
  idx = indexin(labels(m, keys(m, :block)), c)

  x0 = unlist(m)

  iter, p, chains = size(c.data)
  values = Array(Float64, iter, 1, chains)
  for k in 1:chains
    print("\nPROCESSING MCMCChains $(k)/$(chains)\n")
    pct = 0
    for i in 1:iter
      if floor(100 * i / iter) >= pct
        print(string("Row: ", lpad(i, length(string(iter)), ' '),
          "/$(iter) [", lpad(pct, 3, ' '), "%] @ $(strftime(time()))\n"))
        pct += 10
      end
      relist!(m, c.data[i,idx,k][:])
      values[i,1,k] = mapreduce(key -> logpdf(m[key]), +, nkeys)
    end
  end
  print("\n")

  relist!(m, x0)

  values
end

function quantile(c::MCMCChains; q::Vector=[0.025, 0.25, 0.5, 0.75, 0.975])
  X = combine(c)
  dim = size(X)
  labels = map(x -> string(100 * x) * "%", q)
  vals = mapreduce(i -> quantile(X[:,i], q), hcat, 1:dim[2])'
  ChainSummary(vals, c.names, labels, header(c))
end

function summarystats(c::MCMCChains; batchsize::Integer=100)
  X = combine(c)
  n, p = size(X)
  f = (x)->[mean(x), sqrt(var(x)), sem(x), batchSE(x, size=batchsize)]
  labels =["Mean", "SD", "Naive SE", "Batch SE", "ESS"]
  vals = mapreduce(i -> f(X[:,i]), hcat, 1:p)'
  vals = [vals (n * min(vals[:,3] ./ vals[:,4], 1))]
  ChainSummary(vals, c.names, labels, header(c))
end


