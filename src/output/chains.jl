#################### Chains Constructor ####################

function Chains{T<:String}(iters::Integer, params::Integer;
           start::Integer=1, thin::Integer=1, chains::Integer=1,
           names::Vector{T}=String[], model::Model=Model())
  value = Array(Float64, length(start:thin:iters), params, chains)
  fill!(value, NaN)
  Chains(value, start=start, thin=thin, names=names, model=model)
end

function Chains{T<:Real,U<:String,V<:Integer}(value::Array{T,3};
           start::Integer=1, thin::Integer=1, names::Vector{U}=String[],
           chains::Vector{V}=Integer[], model::Model=Model())
  n, p, m = size(value)

  if length(names) == 0
    names = String[string("Param", i) for i in 1:p]
  elseif length(names) != p
    error("size(value, 2) and names dimensions must match")
  end

  if length(chains) == 0
    chains = Integer[1:m;]
  elseif length(chains) != m
    error("size(value, 3) and chains dimensions must match")
  end

  v = convert(Array{Float64, 3}, value)
  Chains(v, range(start, thin, n), String[names...], Integer[chains...], model)
end

function Chains{T<:Real,U<:String}(value::Matrix{T};
           start::Integer=1, thin::Integer=1, names::Vector{U}=String[],
           chains::Integer=1, model::Model=Model())
  Chains(reshape(value, size(value, 1), size(value, 2), 1), start=start,
         thin=thin, names=names, chains=Integer[chains], model=model)
end

function Chains{T<:Real}(value::Vector{T};
           start::Integer=1, thin::Integer=1, names::String="Param1",
           chains::Integer=1, model::Model=Model())
  Chains(reshape(value, length(value), 1, 1), start=start, thin=thin,
         names=String[names], chains=Integer[chains], model=model)
end


#################### Chains Indexing ####################

function Base.getindex(c::Chains, window, names, chains)
  inds1 = window2inds(c, window)
  inds2 = names2inds(c, names)
  Chains(c.value[inds1, inds2, chains],
         start = first(c.range) + (first(inds1) - 1) * step(c.range),
         thin = step(inds1) * step(c.range), names = c.names[inds2],
         chains = c.chains[chains], model = c.model)
end

function Base.setindex!(c::Chains, value, iters, names, chains)
  setindex!(c.value, value, iters2inds(c, iters), names2inds(c, names), chains)
end

window2inds(c::Chains, window) =
  error("indexing Chains iterations with $(typeof(window)) is not supported")
window2inds(c::Chains, ::Colon) = window2inds(c, 1:size(c,1))
window2inds(c::Chains, window::Range) = begin
  range = (window - first(c.range)) / step(c.range) + 1.0
  a = max(ceil(Integer, first(range)), 1)
  b = step(window)
  c = min(floor(Integer, last(range)), size(c.value,1))
  a:b:c
end

iters2inds(c::Chains, iters) = iters
iters2inds(c::Chains, ::Colon) = 1:size(c.value,1)
iters2inds(c::Chains, iters::Range) =
  convert(StepRange{Int,Int}, (iters - first(c.range)) / step(c.range) + 1.0)
iters2inds(c::Chains, iter::Real) =
  Int((iter - first(c.range)) / step(c.range) + 1.0)
iters2inds{T<:Real}(c::Chains, iters::Vector{T}) = begin
  shift, scale = first(c.range), step(c.range)
  Int[(i - shift) / scale + 1.0 for i in iters]
end

names2inds(c::Chains, names) = names
names2inds(c::Chains, ::Colon) = 1:size(c.value,2)
names2inds(c::Chains, name::Real) = [name]
names2inds(c::Chains, name::String) = names2inds(c, [name])
names2inds{T<:String}(c::Chains, names::Vector{T}) = indexin(names, c.names)


#################### Chains Base/Utility Methods ####################

function Base.keys(c::Chains)
  c.names
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
  value = Array(Float64, n * m, p)
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
    "Chains = $(join(map(string, c.chains), ","))\n",
    "Samples per chain = $(length(c.range))\n"
  )
end

function ismodelbased(c::Chains)
  c.model.iter > 0
end

function link(c::Chains)
  cc = copy(c.value)
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
