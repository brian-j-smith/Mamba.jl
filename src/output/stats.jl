#################### Posterior Statistics ####################

function autocor(c::MCMCChains; lags::Vector=[1,5,10,50], relative::Bool=true)
  if relative
    lags *= step(c.range)
  elseif any(lags .% step(c.range) .!= 0)
    error("lags do not correspond to thinning interval")
  end
  labels = map(x -> "Lag " * string(x), lags)
  vals = mapslices(x -> autocor(x, lags)', c.value, [1,2])
  ChainSummary(vals, c.names, labels, header(c))
end

function cor(c::MCMCChains)
  ChainSummary(cor(combine(c)), c.names, c.names, header(c))
end

function describe(c::MCMCChains; q::Vector=[0.025, 0.25, 0.5, 0.75, 0.975],
           etype=:bm, args...)
  println(header(c))
  print("Empirical Posterior Estimates:\n")
  showall(summarystats(c; etype=etype, args...), false)
  print("\nQuantiles:\n")
  showall(quantile(c, q=q), false)
end

function dic(c::MCMCChains)
  m = c.model
  nkeys = keys(m, :output)
  idx = indexin(names(m, keys(m, :block)), c)

  x0 = unlist(m)

  xbar = map(i -> mean(c.value[:,i,:]), idx)
  relist!(m, xbar)
  Dhat = -2.0 * mapreduce(key -> logpdf(m[key]), +, nkeys)
  D = -2.0 * logpdf(c, nkeys)
  p = [mean(D) - Dhat, 0.5 * var(D)]

  relist!(m, x0)

  ChainSummary([Dhat .+ 2.0 * p p], ["pD", "pV"],
               ["DIC", "Effective Parameters"], header(c))
end

function hpd{T<:Real}(x::Vector{T}; alpha::Real=0.05)
  n = length(x)
  m = max(1, ceil(alpha * n))

  y = sort(x)
  a = y[1:m]
  b = y[(n - m + 1):n]
  i = sortperm(b - a)[1]

  [a[i], b[i]]
end

function hpd(c::MCMCChains; alpha::Real=0.05)
  cc = combine(c)
  labels = [string(100 * alpha / 2) * "%", string(100 * (1 - alpha / 2)) * "%"]
  vals = mapslices(x -> hpd(x, alpha=alpha), cc, 1)'
  ChainSummary(vals, c.names, labels, header(c))
end

function logpdf{T<:String}(c::MCMCChains, nkeys::Vector{T})
  m = c.model
  idx = indexin(names(m, keys(m, :block)), c)

  x0 = unlist(m)

  iter, p, chains = size(c.value)
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
      relist!(m, c.value[i,idx,k][:])
      values[i,1,k] = mapreduce(key -> logpdf(m[key]), +, nkeys)
    end
  end
  print("\n")

  relist!(m, x0)

  values
end

function quantile(c::MCMCChains; q::Vector=[0.025, 0.25, 0.5, 0.75, 0.975])
  cc = combine(c)
  labels = map(x -> string(100 * x) * "%", q)
  vals = mapslices(x -> quantile(x, q), cc, 1)'
  ChainSummary(vals, c.names, labels, header(c))
end

function summarystats(c::MCMCChains; etype=:bm, args...)
  cc = combine(c)
  f = x -> [mean(x), std(x), sem(x), mcse(x, etype; args...)]
  labels = ["Mean", "SD", "Naive SE", "MCSE", "ESS"]
  vals = mapslices(x -> f(x), cc, 1)'
  vals = [vals  size(cc, 1) * min(vals[:,3] ./ vals[:,4], 1.0)]
  ChainSummary(vals, c.names, labels, header(c))
end
