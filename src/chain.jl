#################### ChainSummary Type ####################

type ChainSummary
  data::Array{Float64,3}
  rownames::Vector{String}
  colnames::Vector{String}
  header::String

  function ChainSummary(data::Array{Float64,3}, rownames::Vector{String},
                        colnames::Vector{String}, header::String)
    dim = size(data)
    length(rownames) == dim[1] ||
      error("length of rownames not equal to number of rows")
    length(colnames) == dim[2] ||
      error("length of colnames not equal to number of columns")
    new(data, rownames, colnames, header)
  end
end

function ChainSummary{T<:String,U<:String}(data::Array{Float64,3},
           rownames::Vector{T}, colnames::Vector{U}, header::String)
  ChainSummary(deepcopy(data), String[rownames...], String[colnames...], header)
end

function ChainSummary{T<:String,U<:String}(data::Matrix{Float64},
           rownames::Vector{T}, colnames::Vector{U}, header::String)
  dim = size(data)
  ChainSummary(reshape(data, dim[1], dim[2], 1), String[rownames...],
               String[colnames...], header)
end

function annotate(x::Matrix, rownames::Vector, colnames::Vector)
  hcat(["", rownames], vcat(colnames', x))
end

function Base.show(io::IO, s::ChainSummary)
  if size(s.data)[3] == 1
    x = annotate(s.data[:,:,1], s.rownames, s.colnames)
  else
    x = mapslices(x -> annotate(x, s.rownames, s.colnames), s.data, [1,2])
  end
  showall(io, x)
  print("\n")
end

function Base.showall(io::IO, s::ChainSummary)
  println(io, s.header)
  show(io, s)
end


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

function Base.ndims(c::MCMCChain)
  ndims(c.data)
end

function Base.show(io::IO, c::MCMCChain)
  print(io, "Object of type \"$(summary(c))\"\n\n")
  println(io, header(c))
  show(io, c.data)
  print(io, "\n")
end

function Base.size(c::MCMCChain)
  dim = size(c.data)
  (c.start + (dim[1] - 1) * c.thin, dim[2], dim[3])
end

function Base.size(c::MCMCChain, ind)
  size(c)[ind]
end

function combine(c::MCMCChain)
  mapreduce(i -> c.data[:,:,i], vcat, 1:size(c.data)[3])
end

function header(c::MCMCChain)
  dim = size(c.data)
  n = c.start + (dim[1] - 1) * c.thin
  string(
    "Iterations = $(c.start):$n\n",
    "Thinning interval = $(c.thin)\n",
    "Number of chains = $(dim[3])\n",
    "Samples per chain = $(dim[1])\n"
  )
end


#################### MCMCChain Summary Methods ####################

function autocor(c::MCMCChain; lags::Vector=[1,5,10,50], relative::Bool=true)
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

function cor(c::MCMCChain)
  ChainSummary(cor(combine(c)), c.names, c.names, header(c))
end

function describe(c::MCMCChain; batchsize::Integer=100,
                  q::Vector=[0.025, 0.25, 0.5, 0.75, 0.975])
  println(header(c))
  print("Empirical Posterior Estimates:\n")
  show(summarystats(c, batchsize=batchsize))
  print("\nQuantiles:\n")
  show(quantile(c, q=q))
end

function dic(c::MCMCChain)
  m = c.model

  keys = blockkeys(m)
  idx = indexin(labels(m, keys), c)

  x0 = unlist(m, keys)

  xbar = map(i -> mean(c.data[:,i,:]), idx)
  relist!(m, xbar)
  D = -2 * logpdf(m)
  Dbar = -2 * mean(logpdf(c))
  pD = Dbar - D

  relist!(m, x0)

  ChainSummary([D + 2 * pD pD], ["Model"], ["DIC", "pD"], header(c))
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

function hpd(c::MCMCChain; alpha::Real=0.05)
  X = combine(c)
  labels = [string(100 * alpha / 2) * "%", string(100 * (1 - alpha / 2)) * "%"]
  vals = mapreduce(i -> hpd(X[:,i], alpha=alpha), hcat, 1:size(X)[2])
  ChainSummary(vals', c.names, labels, header(c))
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
      values[i,1,k] = logpdf!(m, c.data[i,idx,k][:])
    end
  end
  print("\n")

  relist!(m, x0)

  values
end

function quantile(c::MCMCChain; q::Vector=[0.025, 0.25, 0.5, 0.75, 0.975])
  X = combine(c)
  dim = size(X)
  labels = map(x -> string(100 * x) * "%", q)
  vals = mapreduce(i -> quantile(X[:,i], q), hcat, 1:dim[2])'
  ChainSummary(vals, c.names, labels, header(c))
end

function summarystats(c::MCMCChain; batchsize::Integer=100)
  X = combine(c)
  n, p = size(X)
  f = (x)->[mean(x), sqrt(var(x)), sem(x), batchSE(x, size=batchsize)]
  labels =["Mean", "SD", "Naive SE", "Batch SE", "ESS"]
  vals = mapreduce(i -> f(X[:,i]), hcat, 1:p)'
  vals = [vals (n * min(vals[:,3] ./ vals[:,4], 1))]
  ChainSummary(vals, c.names, labels, header(c))
end


#################### MCMCChain Diagnostic Methods ####################

function gelmandiag(c::MCMCChain; alpha::Real=0.05, mpsrf::Bool=false,
                    transform::Bool=false)
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
  psrf_names = c.names

  if mpsrf
    x = [(n - 1) / n + (m + 1) / m * eigmax(inv(cholfact(W)) * B) / n NaN]
    psrf = vcat(psrf, x)
    psrf_names = [psrf_names, "Multivariate"]
  end

  ChainSummary(psrf, psrf_names, psrf_labels, header(c))
end

function transform_support(c::MCMCChain)
  x = deepcopy(c.data)
  for i in 1:size(x)[2]
    a = minimum(x[:,i,:])
    b = maximum(x[:,i,:])
    if a > 0
      if b < 1
        x[:,i,:] = logit(x[:,i,:])
      else
        x[:,i,:] = log(x[:,i,:])
      end
    end
  end
  x
end
