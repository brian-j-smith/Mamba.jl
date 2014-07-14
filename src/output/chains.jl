#################### Chains Constructor ####################

function Chains{T<:String}(iters::Integer, params::Integer;
           start::Integer=1, thin::Integer=1, chains::Integer=1,
           names::Vector{T}=String[], model::Model=Model())
  value = Array(VariateType, length(start:thin:iters), params, chains)
  fill!(value, NaN)
  Chains(value, start=start, thin=thin, names=names, model=model)
end

function Chains{T<:Real,U<:String}(value::Array{T,3};
           start::Integer=1, thin::Integer=1, names::Vector{U}=String[],
           model::Model=Model())
  n, p, m = size(value)
  if length(names) == 0
    names = String[string("Param", i) for i in 1:p]
  elseif length(names) != p
    error("value column and names dimensions must match")
  end
  vvalue = convert(Array{VariateType, 3}, value)
  Chains(vvalue, String[names...], range(start, thin, n), model)
end

function Chains{T<:Real,U<:String}(value::Matrix{T};
           start::Integer=1, thin::Integer=1, names::Vector{U}=String[],
           model::Model=Model())
  Chains(reshape(value, size(value, 1), size(value, 2), 1), start=start,
         thin=thin, names=names, model=model)
end

function Chains{T<:Real}(value::Vector{T};
           start::Integer=1, thin::Integer=1, names::String="Param1",
           model::Model=Model())
  Chains(reshape(value, length(value), 1, 1), start=start, thin=thin,
         names=String[names], model=model)
end


#################### Chains Base/Utility Methods ####################

function Base.cat(c::Vector{Chains})
  n = length(c)
  n > 0 || error("no Chains to concatenate")
  value = cat(3, [c[i].value for i in 1:n]...)
  range = c[1].range
  Chains(value, start=start(range), thin=step(range), names=c[1].names,
         model=c[1].model)
end

function Base.getindex(c::Chains, iters::Range, names, chains)
  from = max(iceil((first(iters) - first(c.range)) / step(c.range) + 1), 1)
  thin = step(iters)
  to = min(ifloor((last(iters) - first(c.range)) / step(c.range) + 1),
           length(c.range))

  Chains(c.value[from:thin:to, names, chains],
         start = first(c.range) + (from - 1) * step(c.range),
         thin = thin * step(c.range), names = c.names[names],
         model = c.model)
end

function Base.getindex(c::Chains, iters::Range, names::String, chains)
  getindex(c, iters, [names], chains)
end

function Base.getindex{T<:String}(c::Chains, iters::Range, names::Vector{T},
           chains)
  idx = findin(c.names, names)
  length(idx) == length(names) || throw(BoundsError())
  getindex(c, iters, idx, chains)
end

function Base.indexin(names::Vector{String}, c::Chains)
  idx = indexin(names, c.names)
  in(0, idx) && error("node name matches not found in Chains")
  idx
end

function Base.keys(c::Chains)
  c.names
end

function Base.ndims(c::Chains)
  ndims(c.value)
end

function Base.setindex!(c::Chains, value, iters::Range, names, chains)
  setindex!(c, value, collect(iters), names, chains)
end

function Base.setindex!(c::Chains, value, iters::Real, names, chains)
  setindex!(c, value, [iters], names, chains)
end

function Base.setindex!{T<:Real}(c::Chains, value, iters::Vector{T}, names,
           chains)
  start, thin = first(c.range), step(c.range)
  idx = Int[(x - start) / thin + 1.0 for x in iters]
  setindex!(c.value, value, idx, names, chains)
end

function Base.setindex!{T<:Real}(c::Chains, value, iters::Vector{T},
           names::String, chains)
  setindex!(c, value, iters, [names], chains)
end

function Base.setindex!{T<:Real,U<:String}(c::Chains, value, iters::Vector{T},
           names::Vector{U}, chains)
  idx = findin(c.names, names)
  length(idx) == length(names) || throw(BoundsError())
  setindex!(c, value, iters, idx, chains)
end

function Base.show(io::IO, c::Chains)
  print(io, "Object of type \"$(summary(c))\"\n\n")
  println(io, header(c))
  show(io, c.value)
end

function Base.size(c::Chains)
  dim = size(c.value)
  last(c.range), dim[2], dim[3]
end

function Base.size(c::Chains, ind)
  size(c)[ind]
end

function combine(c::Chains)
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

function header(c::Chains)
  string(
    "Iterations = $(first(c.range)):$(last(c.range))\n",
    "Thinning interval = $(step(c.range))\n",
    "Number of chains = $(size(c, 3))\n",
    "Samples per chain = $(length(c.range))\n"
  )
end

function ismodelbased(c::Chains)
  c.model.iter > 0
end

function link(c::Chains)
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
