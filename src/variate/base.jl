#################### Variate Base Methods ####################

for op in [:(Base.endof), :(Base.length), :(Base.ndims), :(Base.size)]
  @eval ($op)(v::Variate) = ($op)(v.data)
end

function Base.convert{T<:Real,N}(::Type{Array{T,N}}, v::VariateVecOrMat)
  convert(Array{T,N}, v.data)
end

function Base.convert{T<:Real}(::Type{T}, v::VariateScalar)
  convert(T, v.data)
end

function Base.getindex(v::VariateVecOrMat, inds...)
  getindex(v.data, inds...)
end

function Base.getindex(v::VariateScalar, i::Real)
  v.data[i]
end

function Base.getindex(v::VariateScalar, inds)
  map(i -> v.data[i], inds)
end

function Base.setindex!(v::VariateVecOrMat, x, inds...)
  setindex!(v.data, x, inds...)
end

function Base.setindex!(v::VariateScalar, x, inds)
  length(x) == 1 || throw(ErrorException("argument dimensions must match"))
  (length(inds) == 1 && collect(inds)[1] == 1) || throw(BoundsError())
  v.data = x[1]
end

function Base.show(io::IO, v::Variate)
  print(io, "Object of type \"$(summary(v))\"\n")
  show(io, v.data)
  print(io, "\n")
end

function Base.showcompact(io::IO, v::Variate)
  showcompact(io, v.data)
end

function names(v::VariateScalar, prefix::String)
  String[prefix]
end

function names(v::VariateVector, prefix::String)
  values = String[]
  for i in 1:length(v)
    push!(values, string(prefix, "[", i, "]"))
  end
  values
end

function names(v::VariateMatrix, prefix::String)
  values = String[]
  for j in 1:size(v, 2)
    for i in 1:size(v, 1)
      push!(values, string(prefix, "[", i, ",", j, "]"))
    end
  end
  values
end
