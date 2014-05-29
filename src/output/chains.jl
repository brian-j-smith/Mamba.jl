#################### MCMCChains Constructor ####################

function MCMCChains{T<:Real,U<:String}(value::Array{T,2}, names::Vector{U};
           start::Integer=1, thin::Integer=1, model::MCMCModel=MCMCModel())
  MCMCChains(reshape(value, size(value, 1), size(value, 2), 1), names,
             start=start, thin=thin, model=model)
end

function MCMCChains{T<:Real,U<:String}(value::Array{T,3}, names::Vector{U};
           start::Integer=1, thin::Integer=1, model::MCMCModel=MCMCModel())
  length(names) == size(value, 2) ||
    error("value column and names dimensions must match")
  varval = convert(Array{VariateType, 3}, value)
  MCMCChains(varval, String[names...], range(start, thin, size(varval, 1)),
             model)
end

function MCMCChains{T<:String}(iter::Integer, names::Vector{T};
           start::Integer=1, thin::Integer=1, chains::Integer=1,
           model::MCMCModel=MCMCModel())
  varval = Array(VariateType, iter, length(names), chains)
  fill!(varval, NaN)
  MCMCChains(varval, String[names...], start=start, thin=thin, model=model)
end


#################### MCMCChains Base/Utility Methods ####################

function Base.getindex{T<:String}(c::MCMCChains, iter::Range, names::Vector{T},
           chains::Vector)
  from = max(iceil((first(iter) - first(c.range)) / step(c.range) + 1), 1)
  thin = step(iter)
  to = min(ifloor((last(iter) - first(c.range)) / step(c.range) + 1),
           length(c.range))

  idx1 = from:thin:to
  idx2 = findin(c.names, names)
  idx3 = findin(1:size(c, 3), chains)

  value = c.value[idx1, idx2, idx3]
  MCMCChains(value, c.names[idx2],
             range(first(c.range) + (from - 1) * step(c.range),
                   thin * step(c.range), length(idx1)),
             c.model)
end

function Base.getindex(c::MCMCChains, iter::Range, names::Vector, chains::Vector)
  c[iter, c.names[names], chains]
end

function Base.getindex(c::MCMCChains, inds...)
  length(inds) == 3 ||
    error("must supply 3-dimensional index for iter, names, and chains")
  idx = inds[1]
  iter = isa(idx, Range1) ? (first(idx):1:last(idx)) :
         isa(idx, Range) ? idx : throw(TypeError())
  names = collect(inds[2])
  chains = collect(inds[3])
  c[iter, names, chains]
end

function Base.indexin(names::Vector{String}, c::MCMCChains)
  idx = indexin(names, c.names)
  all(idx .!= 0) || error("node name matches not found in MCMCChains")
  idx
end

function Base.keys(c::MCMCChains)
  c.names
end

function Base.ndims(c::MCMCChains)
  ndims(c.value)
end

function Base.show(io::IO, c::MCMCChains)
  print(io, "Object of type \"$(summary(c))\"\n\n")
  println(io, header(c))
  show(io, c.value)
  print(io, "\n")
end

function Base.size(c::MCMCChains)
  dim = size(c.value)
  last(c.range), dim[2], dim[3]
end

function Base.size(c::MCMCChains, ind)
  size(c)[ind]
end

function combine(c::MCMCChains)
  n, p, m = size(c.value)
  value = Array(VariateType, n * m, p)
  for j in 1:p
    idx = 1
    for i in 1:n, k in 1:m
      value[idx, j] = c.value[i, j, k]
      idx += 1
    end
  end
  value
end

function header(c::MCMCChains)
  string(
    "Iterations = $(first(c.range)):$(last(c.range))\n",
    "Thinning interval = $(step(c.range))\n",
    "Number of chains = $(size(c, 3))\n",
    "Samples per chain = $(length(c.range))\n"
  )
end

function link(c::MCMCChains)
  cc = deepcopy(c.value)
  idx0 = 1:length(c.names)
  for key in intersect(keys(c.model, :monitor), keys(c.model, :stochastic))
    node = c.model[key]
    idx = nonzeros(indexin(names(node), c.names))
    if length(idx) > 0
      cc[:,idx,:] = mapslices(x -> link(node, x), cc[:,idx,:], 2)
      idx0 = setdiff(idx0, idx)
    end
  end
  for j in idx0
    x = cc[:,j,:]
    if minimum(x) > 0.0
      cc[:,j,:] = maximum(x) < 1.0 ? logit(x) : log(x)
    end
  end
  cc
end
