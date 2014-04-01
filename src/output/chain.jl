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
