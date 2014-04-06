#################### MCMCChains Constructor ####################

function MCMCChains{T<:String}(names::Vector{T}, iter::Integer;
           start::Integer=1, thin::Integer=1, chains::Integer=1,
           model::MCMCModel=MCMCModel())
  data = Array(VariateType, iter, length(names), chains)
  fill!(data, NaN)
  MCMCChains(data, String[names...], start, thin, model)
end


#################### MCMCChains Base/Utility Methods ####################

function Base.getindex{T<:String}(c::MCMCChains, iter::Range, names::Vector{T},
           chains::Vector)
  dim = size(c.data)

  from = max(iceil((first(iter) - c.start) / c.thin + 1), 1)
  thin = step(iter)
  to = min(ifloor((last(iter) - c.start) / c.thin + 1), dim[1])

  idx1 = from:thin:to
  idx2 = findin(c.names, names)
  idx3 = findin(1:dim[3], chains)

  data = c.data[idx1, idx2, idx3]
  MCMCChains(data, c.names[idx2], c.start + (from - 1) * c.thin, c.thin * thin,
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
  ndims(c.data)
end

function Base.show(io::IO, c::MCMCChains)
  print(io, "Object of type \"$(summary(c))\"\n\n")
  println(io, header(c))
  show(io, c.data)
  print(io, "\n")
end

function Base.size(c::MCMCChains)
  dim = size(c.data)
  (c.start + (dim[1] - 1) * c.thin, dim[2], dim[3])
end

function Base.size(c::MCMCChains, ind)
  size(c)[ind]
end

function combine(c::MCMCChains)
  mapreduce(i -> c.data[:,:,i], vcat, 1:size(c.data, 3))
end

function header(c::MCMCChains)
  dim = size(c.data)
  n = c.start + (dim[1] - 1) * c.thin
  string(
    "Iterations = $(c.start):$n\n",
    "Thinning interval = $(c.thin)\n",
    "Number of chains = $(dim[3])\n",
    "Samples per chain = $(dim[1])\n"
  )
end

function link(c::MCMCChains)
  X = deepcopy(c.data)
  m = size(X, 1)
  for key in keys(c.model, true)
    node = c.model[key]
    idx = nonzeros(indexin(labels(c.model, [key]), c.names))
    if length(idx) > 0
      if isa(node, MCMCStochastic)
        X[:,idx,:] = mapslices(x -> link(node.distr, x), X[:,idx,:], 2)
      else
        for j in idx
          x = X[:,j,:]
          if minimum(x) > 0.0
            X[:,j,:] = maximum(x) < 1.0 ? logit(x) : log(x)
          end
        end
      end
    end
  end
  X
end
