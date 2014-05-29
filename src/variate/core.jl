#################### Variate Core Methods ####################

for op in [:(Base.endof), :(Base.length), :(Base.ndims), :(Base.size)]
  @eval ($op)(v::AbstractVariate) = ($op)(v.value)
end

Base.size(v::AbstractVariate, d) = Base.size(v)[d]

function Base.convert{T<:Real,N}(::Type{Array{T,N}}, v::AbstractVariateArray{N})
  convert(Array{T,N}, v.value)
end

function Base.convert{T<:Real}(::Type{T}, v::AbstractVariateScalar)
  convert(T, v.value)
end

function Base.getindex{N}(v::AbstractVariateArray{N}, inds...)
  getindex(v.value, inds...)
end

function Base.getindex(v::AbstractVariateScalar, i::Real)
  v.value[i]
end

function Base.getindex(v::AbstractVariateScalar, inds)
  map(i -> v.value[i], inds)
end

function Base.setindex!{N}(v::AbstractVariateArray{N}, x, inds...)
  setindex!(v.value, x, inds...)
end

function Base.setindex!(v::AbstractVariateScalar, x, inds)
  length(x) == 1 || error("value to store must be of length 1")
  collect(Int, inds) == [1] || throw(BoundsError())
  v.value = x[]
end

function Base.show(io::IO, v::AbstractVariate)
  print(io, "Object of type \"$(summary(v))\"\n")
  show(io, v.value)
  print(io, "\n")
end

function Base.showcompact(io::IO, v::AbstractVariate)
  showcompact(io, v.value)
end

function names(v::AbstractVariateScalar, prefix)
  String[string(prefix)]
end

function names{N}(v::AbstractVariateArray{N}, prefix)
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
