#################### Posterior Statistics ####################

function autocor(c::Chains; lags::Vector=[1,5,10,50], relative::Bool=true)
  if relative
    lags *= step(c.range)
  elseif any(lags .% step(c.range) .!= 0)
    error("lags do not correspond to thinning interval")
  end
  labels = map(x -> "Lag " * string(x), lags)
  vals = mapslices(x -> autocor(x, lags)', c.value, [1,2])
  ChainSummary(vals, c.names, labels, header(c))
end

function cor(c::Chains)
  ChainSummary(cor(combine(c)), c.names, c.names, header(c))
end

describe(c::Chains; args...) = describe(STDOUT, c; args...)

function describe(io::IO, c::Chains; q::Vector=[0.025, 0.25, 0.5, 0.75, 0.975],
                  etype=:bm, args...)
  ps_stats = summarystats(c; etype=etype, args...)
  ps_quantiles = quantile(c, q=q)
  println(io, ps_stats.header)
  print(io, "Empirical Posterior Estimates:\n")
  show(io, ps_stats)
  print(io, "Quantiles:\n")
  show(io, ps_quantiles)
end

function dic(c::Chains)
  ismodelbased(c) || error("dic requires Chains from a Model fit")

  m = c.model
  nkeys = keys(m, :output)
  idx = indexin(names(m, keys(m, :block)), c.names)
  in(0, idx) && error("dic requires all sampled nodes to be monitored")

  xbar = map(i -> mean(c.value[:,i,:]), idx)
  relist!(m, xbar)
  Dhat = -2.0 * mapreduce(key -> logpdf(m[key]), +, nkeys)
  D = -2.0 * logpdf(c, nkeys)
  p = [mean(D) - Dhat, 0.5 * var(D)]

  ChainSummary([Dhat + 2.0 * p  p], ["pD", "pV"],
               ["DIC", "Effective Parameters"], header(c))
end

function hpd{T<:Real}(x::Vector{T}; alpha::Real=0.05)
  n = length(x)
  m = max(1, ceil(Integer, alpha * n))

  y = sort(x)
  a = y[1:m]
  b = y[(n - m + 1):n]
  _,i = findmin(b - a)

  [a[i], b[i]]
end

function hpd(c::Chains; alpha::Real=0.05)
  cc = combine(c)
  pct = first(showoff([100.0 * (1.0 - alpha)]))
  labels = ["$(pct)% Lower", "$(pct)% Upper"]
  vals = mapslices(x -> hpd(x, alpha=alpha), cc, 1)'
  ChainSummary(vals, c.names, labels, header(c))
end

function logpdf(c::Chains, nkeys::Vector{Symbol})
  ismodelbased(c) || error("logpdf requires Chains from a Model fit")

  m = c.model
  idx = indexin(names(m, keys(m, :block)), c.names)
  in(0, idx) && error("logpdf requires all sampled nodes to be monitored")

  iters, p, chains = size(c.value)
  values = Array(Float64, iters, 1, chains)
  frame = ChainProgressFrame(
    "MCMC Processing of $iters Iterations x $chains Chain" * "s"^(chains > 1),
    true
  )
  for k in 1:chains
    meter = ChainProgress(frame, k, iters)
    for i in 1:iters
      relist!(m, vec(c.value[i,idx,k]))
      values[i,1,k] = mapreduce(key -> logpdf(m[key]), +, nkeys)
      next!(meter)
    end
    println()
  end

  values
end

function quantile(c::Chains; q::Vector=[0.025, 0.25, 0.5, 0.75, 0.975])
  cc = combine(c)
  labels = map(x -> string(100 * x) * "%", q)
  vals = mapslices(x -> quantile(x, q), cc, 1)'
  ChainSummary(vals, c.names, labels, header(c))
end

function summarystats(c::Chains; etype=:bm, args...)
  cc = combine(c)
  f = x -> [mean(x), std(x), sem(x), mcse(x, etype; args...)]
  labels = ["Mean", "SD", "Naive SE", "MCSE", "ESS"]
  vals = mapslices(x -> f(x), cc, 1)'
  vals = [vals  min((vals[:,2] ./ vals[:,4]).^2, size(cc, 1))]
  ChainSummary(vals, c.names, labels, header(c))
end
