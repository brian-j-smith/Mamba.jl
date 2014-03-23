#################### MCMCChain Constructor ####################

function MCMCChain(params::Vector, iter::Integer; start::Integer=1,
                   thin::Integer=1, chains::Integer=1,
                   model::MCMCModel=MCMCModel())
  data = Array(VariateType, iter, length(params), chains)
  fill!(data, NaN)
  MCMCChain(data, String[params...], start, thin, model)
end


#################### MCMCChain Base/Utility Methods ####################

function Base.getindex{T<:String}(c::MCMCChain, iter::Range, names::Vector{T},
                                  chains::Vector)
  dim = size(c.data)

  from = max(iceil((first(iter) - c.start) / c.thin + 1), 1)
  thin = step(iter)
  to = min(ifloor((last(iter) - c.start) / c.thin + 1), dim[1])

  idx1 = from:thin:to
  idx2 = findin(c.names, names)
  idx3 = findin(1:dim[3], chains)

  data = c.data[idx1, idx2, idx3]
  MCMCChain(data, c.names[idx2], c.start + (from - 1) * c.thin, c.thin * thin,
            c.model)
end

function Base.getindex(c::MCMCChain, iter::Range, names::Vector, chains::Vector)
  c[iter, c.names[names], chains]
end

function Base.getindex(c::MCMCChain, inds...)
  length(inds) == 3 ||
    error("must supply 3-dimensional index for iter, names, and chains")
  idx = inds[1]
  iter = isa(idx, Range1) ? (first(idx):1:last(idx)) :
         isa(idx, Range) ? idx : throw(TypeError())
  names = collect(inds[2])
  chains = collect(inds[3])
  c[iter, names, chains]
end

function Base.indexin(names::Vector{String}, c::MCMCChain)
  idx = indexin(names, c.names)
  all(idx .!= 0) || error("node name matches not found in MCMCChain")
  idx
end

function Base.keys(c::MCMCChain)
  c.names
end

function Base.map(f::Function, c::MCMCChain)
  x = []
  iter, p, chains = size(c.data)
  for k in 1:chains
    print("\nPROCESSING MCMCChain $(k)/$(chains)\n")
    global pct = 0
    g = function(i)
      if floor(100 * i / iter) >= pct
        print(string("Row: ", lpad(i, length(string(iter)), ' '),
          "/$(iter) [", lpad(pct, 3, ' '), "%] @ $(strftime(time()))\n"))
        pct += 10
      end
      f(c.data[i,:,k][:])
    end
    x = [x, map(g, 1:iter)]
  end
  print("\n")
  x
end

function Base.ndims(c::MCMCChain)
  ndims(c.data)
end

function Base.show(io::IO, c::MCMCChain)
  dim = size(c.data)
  n = c.start + (dim[1] - 1) * c.thin
  print(io, "Object of type \"$(summary(c))\"\n\n")
  print(io, string(
    "Iterations = $(c.start):$n\n",
    "Thinning interval = $(c.thin)\n",
    "Number of chains = $(dim[3])\n",
    "Samples per chain = $(dim[1])\n"
  ))
  stats, quants = describe(c)
  print(io, "\nEmpirical Posterior Estimates:\n")
  show(io, stats)
  print(io, "\n\nQuantiles:\n")
  show(io, quants)
  print(io, "\n")
end

function Base.size(c::MCMCChain)
  dim = size(c.data)
  (c.start + (dim[1] - 1) * c.thin, dim[2], dim[3])
end

function Base.size(c::MCMCChain, ind)
  size(c)[ind]
end

function annotate(x::Matrix, colnames::Vector, rownames::Vector, index="")
  hcat([index, rownames], vcat(colnames', x))
end

function combine(c::MCMCChain)
  mapreduce(i -> c.data[:,:,i], vcat, 1:size(c.data)[3])
end


#################### MCMCChain Summary Methods ####################

function Base.cor(c::MCMCChain)
  vals = cor(combine(c))
  annotate(vals, c.names, c.names, "Node")
end

function autocor(c::MCMCChain; lags::Vector=[1,5,10,50], relative::Bool=true)
  if relative
    lags *= c.thin
  elseif any(lags .% c.thin .!= 0)
    error("lags do not correspond to thinning interval")
  end
  labels = map(x -> "Lag " * string(x), lags)
  mapslices(x -> annotate(autocor(x, lags)', labels, c.names, "Node"), c.data, [1,2])
end

function batchSE(x::Vector, size::Integer=100)
  m = div(length(x), size)
  m >= 2 || error("2 or more batches needed to compute SE")
  mbar = [mean(x[i*size+(1:size)]) for i in 0:m-1]
  sem(mbar)
end

function describe(c::MCMCChain; batchsize::Integer=100,
                  q::Vector=[0.025, 0.25, 0.5, 0.75, 0.975])
  X = combine(c)
  n, p = size(X)

  f = (x)->[mean(x), sqrt(var(x)), sem(x), batchSE(x, batchsize)]
  labels =["Mean", "SD", "Naive SE", "Batch SE", "ESS"]
  vals = mapreduce(i -> f(X[:,i]), hcat, 1:p)'
  vals = [vals (n * min(vals[:,3] ./ vals[:,4], 1))]
  stats = annotate(vals, labels, c.names, "Node")

  labels = map(x -> string(100 * x) * "%", q)
  vals = mapreduce(i -> quantile(X[:,i], q), hcat, 1:p)'
  quants = annotate(vals, labels, c.names, "Node")

  stats, quants
end

function dic(c::MCMCChain)
  m = c.model

  keys = blockkeys(m)
  idx = indexin(labels(m, keys), c)

  x0 = unlist(m, keys)

  xbar = map(i -> mean(c.data[:,i,:]), idx)
  relist!(m, xbar, keys)
  update!(m)
  D = -2 * logpdf(m)
  Dbar = -2 * mean(logpdf(c))
  pD = Dbar - D

  relist!(m, x0, keys)
  update!(m)

  annotate([D + 2 * pD pD], ["DIC", "pD"], ["Model"])
end

function hpd(x::Vector, alpha::Real=0.05)
  n = length(x)
  m = max(1, ceil(alpha * n))

  y = sort(x)
  a = y[1:m]
  b = y[(n - m + 1):n]
  i = sortperm(b - a)[1]

  [a[i], b[i]]
end

function hpd(c::MCMCChain, alpha::Real=0.05)
  X = combine(c)
  labels = [string(100 * alpha / 2) * "%", string(100 * (1 - alpha / 2)) * "%"]
  vals = mapreduce(i -> hpd(X[:,i], alpha), hcat, 1:size(X)[2])
  annotate(vals', labels, c.names, "Node")
end

function logpdf(c::MCMCChain)
  m = c.model

  keys = blockkeys(m)
  idx = indexin(labels(m, keys), c)

  x0 = unlist(m, keys)

  iter, p, chains = size(c.data)
  values = Array(Float64, iter, 1, chains)
  for k in 1:chains
    print("\nPROCESSING MCMCChain $(k)/$(chains)\n")
    pct = 0
    for i in 1:iter
      if floor(100 * i / iter) >= pct
        print(string("Row: ", lpad(i, length(string(iter)), ' '),
          "/$(iter) [", lpad(pct, 3, ' '), "%] @ $(strftime(time()))\n"))
        pct += 10
      end
      x = c.data[i,idx,k][:]
      relist!(m, x, keys)
      update!(m)
      values[i,1,k] = logpdf(m)
    end
  end
  print("\n")

  relist!(m, x0, keys)
  update!(m)

  values
end


#################### MCMCChain Diagnostic Methods ####################

function gelmandiag(c::MCMCChain; alpha::Real=0.05, transform::Bool=false)
  n, p, m = size(c.data)
  m >= 2 || error("2 or more chains needed to run gelman diagnostic")

  psi = transform ? transform_support(c) : c.data

  S2 = mapslices(x -> reshape(crosscov(x, x, [0]), p, p), psi, [1,2])
  W = mapreduce(i -> S2[:,:,i], +, 1:m) / m

  psibar = reshape(mapslices(mean, psi, 1), p, m)'
  B = n * reshape(crosscov(psibar, psibar, [0]), p, p)

  w = diag(W)
  b = diag(B)
  s2 = reshape(mapslices(diag, S2, [1,2]), p, m)'
  psibar2 = map(i -> mean(psibar[:,i]), 1:p)

  var_w = map(i -> var(s2[:,i]), 1:p) / m
  var_b = (2 / (m - 1)) * b.^2
  var_wb = (n / m) * (diag(reshape(crosscov(s2, psibar.^2, [0]), p, p)) -
           2 * psibar2 .* diag(reshape(crosscov(s2, psibar, [0]), p, p)))

  V = ((n - 1) / n) * w + ((m + 1) / (m * n)) * b
  var_V = ((n - 1)^2 * var_w + ((m + 1) / m)^2 * var_b +
           (2 * (n - 1) * (m + 1) / m) * var_wb) / n^2
  df = 2 * V.^2 ./ var_V
  B_df = m - 1
  W_df = 2 * w.^2 ./ var_w

  R_fixed = (n - 1) / n
  R_random = ((m + 1) / (m * n)) * b ./ w
  R_est = R_fixed + R_random
  q = 1 - alpha / 2
  R_upper = R_fixed + R_random .*
            map(df2 -> quantile(FDist(B_df, df2), q), W_df)

  psrf = sqrt((df + 3) / (df + 1) * [R_est R_upper])
  psrf_labels = ["PSRF", string(100 * q) * "%"]

  mpsrf = (n - 1) / n + (m + 1) / m * eigmax(inv(cholfact(W)) * B) / n

  annotate(psrf, psrf_labels, c.names, "Node"),
  ["" "MPSRF"; "All Nodes" mpsrf]
end

function transform_support(c::MCMCChain)
  x = deepcopy(c.data)
  for i in 1:size(x)[2]
    a = minimum(x[:,i,:])
    b = maximum(x[:,i,:])
    if a > 0
      if b < 1
        x[:,i,:] = -log(1 / x[:,i,:] - 1)
      else
        x[:,i,:] = log(x[:,i,:])
      end
    end
  end
  x
end
