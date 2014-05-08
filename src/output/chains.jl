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
  MCMCChains(varval, String[names...], start, thin, model)
end

function MCMCChains{T<:String}(iter::Integer, names::Vector{T};
           start::Integer=1, thin::Integer=1, chains::Integer=1,
           model::MCMCModel=MCMCModel())
  varval = Array(VariateType, iter, length(names), chains)
  fill!(varval, NaN)
  MCMCChains(varval, String[names...], start, thin, model)
end


#################### MCMCChains Base/Utility Methods ####################

function Base.getindex{T<:String}(c::MCMCChains, iter::Range, names::Vector{T},
           chains::Vector)
  dim = size(c.value)

  from = max(iceil((first(iter) - c.start) / c.thin + 1), 1)
  thin = step(iter)
  to = min(ifloor((last(iter) - c.start) / c.thin + 1), dim[1])

  idx1 = from:thin:to
  idx2 = findin(c.names, names)
  idx3 = findin(1:dim[3], chains)

  value = c.value[idx1, idx2, idx3]
  MCMCChains(value, c.names[idx2], c.start + (from - 1) * c.thin, c.thin * thin,
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
  (c.start + (dim[1] - 1) * c.thin, dim[2], dim[3])
end

function Base.size(c::MCMCChains, ind)
  size(c)[ind]
end

function combine(c::MCMCChains)
  mapreduce(i -> c.value[:,:,i], vcat, 1:size(c.value, 3))
end

function header(c::MCMCChains)
  dim = size(c.value)
  n = c.start + (dim[1] - 1) * c.thin
  string(
    "Iterations = $(c.start):$n\n",
    "Thinning interval = $(c.thin)\n",
    "Number of chains = $(dim[3])\n",
    "Samples per chain = $(dim[1])\n"
  )
end

function link(c::MCMCChains)
  X = deepcopy(c.value)
  idx0 = 1:length(c.names)
  for key in intersect(keys(c.model, :monitor), keys(c.model, :stochastic))
    node = c.model[key]
    idx = nonzeros(indexin(names(node), c.names))
    if length(idx) > 0
      X[:,idx,:] = mapslices(x -> link(node, x), X[:,idx,:], 2)
      idx0 = setdiff(idx0, idx)
    end
  end
  for j in idx0
    x = X[:,j,:]
    if minimum(x) > 0.0
      X[:,j,:] = maximum(x) < 1.0 ? logit(x) : log(x)
    end
  end
  X
end
