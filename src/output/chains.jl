#################### MCMCChains Constructor ####################

function MCMCChains{T<:Real,U<:String}(value::Array{T,2};
           start::Integer=1, thin::Integer=1, names::Vector{U}=String[],
           model::MCMCModel=MCMCModel())
  MCMCChains(reshape(value, size(value, 1), size(value, 2), 1), start=start,
             thin=thin, names=names, model=model)
end

function MCMCChains{T<:Real,U<:String}(value::Array{T,3};
           start::Integer=1, thin::Integer=1, names::Vector{U}=String[],
           model::MCMCModel=MCMCModel())
  n, p, m = size(value)
  if length(names) == 0
    names = String[string("Param", i) for i in 1:p]
  elseif length(names) != p
    error("value column and names dimensions must match")
  end
  vvalue = convert(Array{VariateType, 3}, value)
  MCMCChains(vvalue, String[names...], range(start, thin, n), model)
end

function MCMCChains{T<:String}(iters::Integer, params::Integer;
           start::Integer=1, thin::Integer=1, chains::Integer=1,
           names::Vector{T}=String[], model::MCMCModel=MCMCModel())
  value = Array(VariateType, length(start:thin:iters), params, chains)
  fill!(value, NaN)
  MCMCChains(value, start=start, thin=thin, names=names, model=model)
end


#################### MCMCChains Base/Utility Methods ####################

function Base.getindex{T<:String}(c::MCMCChains, iters::Range, names::Vector{T},
           chains::Vector)
  from = max(iceil((first(iters) - first(c.range)) / step(c.range) + 1), 1)
  thin = step(iters)
  to = min(ifloor((last(iters) - first(c.range)) / step(c.range) + 1),
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

function Base.getindex(c::MCMCChains, iters::Range, names::Vector, chains::Vector)
  c[iters, c.names[names], chains]
end

function Base.getindex(c::MCMCChains, inds...)
  length(inds) == 3 ||
    error("must supply 3-dimensional index for iters, names, and chains")
  idx = inds[1]
  iters = isa(idx, Range1) ? (first(idx):1:last(idx)) :
          isa(idx, Range) ? idx : throw(TypeError())
  names = collect(inds[2])
  chains = collect(inds[3])
  c[iters, names, chains]
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
    idx = findin(c.names, names(node))
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
