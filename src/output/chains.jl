#################### Chains Methods ####################

#################### Constructors ####################

function Chains{T<:String}(iters::Integer, params::Integer;
           start::Integer=1, thin::Integer=1, chains::Integer=1,
           names::Vector{T}=String[])
  value = Array(Float64, length(start:thin:iters), params, chains)
  fill!(value, NaN)
  Chains(value, start=start, thin=thin, names=names)
end

function Chains{T<:Real,U<:String,V<:Integer}(value::Array{T,3};
           start::Integer=1, thin::Integer=1, names::Vector{U}=String[],
           chains::Vector{V}=Int[])
  n, p, m = size(value)

  if length(names) == 0
    names = String[string("Param", i) for i in 1:p]
  elseif length(names) != p
    error("size(value, 2) and names dimensions must match")
  end

  if length(chains) == 0
    chains = Int[1:m;]
  elseif length(chains) != m
    error("size(value, 3) and chains dimensions must match")
  end

  v = convert(Array{Float64, 3}, value)
  Chains(v, range(start, thin, n), String[names...], Int[chains...])
end

function Chains{T<:Real,U<:String}(value::Matrix{T};
           start::Integer=1, thin::Integer=1, names::Vector{U}=String[],
           chains::Integer=1)
  Chains(reshape(value, size(value, 1), size(value, 2), 1), start=start,
         thin=thin, names=names, chains=Int[chains])
end

function Chains{T<:Real}(value::Vector{T};
           start::Integer=1, thin::Integer=1, names::String="Param1",
           chains::Integer=1)
  Chains(reshape(value, length(value), 1, 1), start=start, thin=thin,
         names=String[names], chains=Int[chains])
end


#################### Indexing ####################

function Base.getindex(c::Chains, window, names, chains)
  inds1 = window2inds(c, window)
  inds2 = names2inds(c, names)
  Chains(c.value[inds1, inds2, chains],
         start = first(c.range) + (first(inds1) - 1) * step(c.range),
         thin = step(inds1) * step(c.range), names = c.names[inds2],
         chains = c.chains[chains])
end

function Base.setindex!(c::AbstractChains, value, iters, names, chains)
  setindex!(c.value, value, iters2inds(c, iters), names2inds(c, names), chains)
end

macro mapiters(iters, c)
  quote
    ($iters - first(($c).range)) / step(($c).range) + 1.0
  end
end

window2inds(c::AbstractChains, window) =
  error("indexing Chains iterations with $(typeof(window)) is not supported")
window2inds(c::AbstractChains, ::Colon) = window2inds(c, 1:size(c,1))
window2inds(c::AbstractChains, window::Range) = begin
  range = @mapiters(window, c)
  a = max(ceil(Int, first(range)), 1)
  b = step(window)
  c = min(floor(Int, last(range)), size(c.value,1))
  a:b:c
end

iters2inds(c::AbstractChains, iters) = iters
iters2inds(c::AbstractChains, ::Colon) = 1:size(c.value,1)
iters2inds(c::AbstractChains, iters::Range) =
  convert(StepRange{Int,Int}, @mapiters(iters, c))
iters2inds(c::AbstractChains, iter::Real) = Int(@mapiters(iter, c))
iters2inds{T<:Real}(c::AbstractChains, iters::Vector{T}) =
  Int[@mapiters(i, c) for i in iters]

names2inds(c::AbstractChains, names) = names
names2inds(c::AbstractChains, ::Colon) = 1:size(c.value,2)
names2inds(c::AbstractChains, name::Real) = [name]
names2inds(c::AbstractChains, name::String) = names2inds(c, [name])
names2inds{T<:String}(c::AbstractChains, names::Vector{T}) =
  indexin(names, c.names)


#################### Auxilliary Functions ####################

function Base.keys(c::AbstractChains)
  c.names
end

function Base.show(io::IO, c::AbstractChains)
  print(io, "Object of type \"$(summary(c))\"\n\n")
  println(io, header(c))
  show(io, c.value)
end

function Base.size(c::AbstractChains)
  dim = size(c.value)
  last(c.range), dim[2], dim[3]
end

function Base.size(c::AbstractChains, ind)
  size(c)[ind]
end

function combine(c::AbstractChains)
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

function header(c::AbstractChains)
  string(
    "Iterations = $(first(c.range)):$(last(c.range))\n",
    "Thinning interval = $(step(c.range))\n",
    "Chains = $(join(map(string, c.chains), ","))\n",
    "Samples per chain = $(length(c.range))\n"
  )
end

function indiscretesupport(c::AbstractChains, bounds::Tuple{Real,Real}=(0,Inf))
  nrows, nvars, nchains = size(c.value)
  result = Array(Bool, nvars * (nrows > 0))
  for i in 1:nvars
    result[i] = true
    for j in 1:nrows, k in 1:nchains
      x = c.value[j,i,k]
      if !isinteger(x) || x < bounds[1] || x > bounds[2]
        result[i] = false
        break
      end
    end
  end
  result
end

function link(c::AbstractChains)
  cc = copy(c.value)
  for j in 1:length(c.names)
    x = cc[:,j,:]
    if minimum(x) > 0.0
      cc[:,j,:] = maximum(x) < 1.0 ? logit(x) : log(x)
    end
  end
  cc
end
