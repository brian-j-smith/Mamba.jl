#################### Variate Core Methods ####################

for op in [:(Base.endof), :(Base.length), :(Base.ndims), :(Base.size)]
  @eval ($op)(v::Variate) = ($op)(v.value)
end

Base.size(v::Variate, d) = Base.size(v)[d]

function Base.convert{T<:Real,N}(::Type{Array{T,N}}, v::MultiVariate{N})
  convert(Array{T,N}, v.value)
end

function Base.convert{T<:Real}(::Type{T}, v::UniVariate)
  convert(T, v.value)
end

function Base.getindex{N}(v::MultiVariate{N}, inds...)
  getindex(v.value, inds...)
end

function Base.getindex(v::UniVariate, i::Real)
  v.value[i]
end

function Base.getindex(v::UniVariate, inds)
  map(i -> v.value[i], inds)
end

function Base.setindex!{N}(v::MultiVariate{N}, x, inds...)
  setindex!(v.value, x, inds...)
end

function Base.setindex!(v::UniVariate, x, inds)
  length(x) == 1 || error("value to store must be of length 1")
  collect(Int, inds) == [1] || throw(BoundsError())
  v.value = x[]
end

function Base.show(io::IO, v::Variate)
  print(io, "Object of type \"$(summary(v))\"\n")
  show(io, v.value)
end

function Base.showcompact(io::IO, v::Variate)
  showcompact(io, v.value)
end

function names(v::UniVariate, prefix)
  String[string(prefix)]
end

function names{N}(v::MultiVariate{N}, prefix)
  n = length(v)
  values = Array(String, n)
  dims = size(v)
  offset = ndims(v) > 1 ? 1 : 2
  for i in 1:n
    s = string(ind2sub(dims, i))
    values[i] = string(prefix, "[", s[2:end-offset], "]")
  end
  values
end
