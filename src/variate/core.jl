#################### Variate Core Methods ####################

for op in [:(Base.endof), :(Base.length), :(Base.ndims), :(Base.size)]
  @eval ($op)(v::Variate) = ($op)(v.value)
end

Base.size(v::Variate, d) = Base.size(v)[d]

function Base.convert{T<:Real,N}(::Type{Array{T,N}}, v::VariateArray{N})
  convert(Array{T,N}, v.value)
end

function Base.convert{T<:Real}(::Type{T}, v::VariateScalar)
  convert(T, v.value)
end

function Base.getindex{N}(v::VariateArray{N}, inds...)
  getindex(v.value, inds...)
end

function Base.getindex(v::VariateScalar, i::Real)
  v.value[i]
end

function Base.getindex(v::VariateScalar, inds)
  map(i -> v.value[i], inds)
end

function Base.setindex!{N}(v::VariateArray{N}, x, inds...)
  setindex!(v.value, x, inds...)
end

function Base.setindex!(v::VariateScalar, x, inds)
  length(x) == 1 || throw(ErrorException("argument dimensions must match"))
  (length(inds) == 1 && collect(inds)[1] == 1) || throw(BoundsError())
  v.value = x[1]
end

function Base.show(io::IO, v::Variate)
  print(io, "Object of type \"$(summary(v))\"\n")
  show(io, v.value)
  print(io, "\n")
end

function Base.showcompact(io::IO, v::Variate)
  showcompact(io, v.value)
end

function names{T<:String}(v::VariateScalar, prefix::T)
  String[prefix]
end

function names{N,T<:String}(v::VariateArray{N}, prefix::T)
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
