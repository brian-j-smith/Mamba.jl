#################### Variate ####################

#################### Conversions ####################

Base.convert(::Type{Bool}, v::ScalarVariate) = convert(Bool, v.value)
Base.convert{T<:Integer}(::Type{T}, v::ScalarVariate) = convert(T, v.value)
Base.convert{T<:AbstractFloat}(::Type{T}, v::ScalarVariate) = convert(T, v.value)

Base.convert{T<:Real,N}(::Union{Type{AbstractArray{T,N}},Type{Array{T,N}}},
                        v::ArrayVariate{N}) = convert(Array{T,N}, v.value)

Base.unsafe_convert{T<:Real}(::Type{Ptr{T}}, v::ArrayVariate) = pointer(v.value)


#################### Base Functions ####################

Base.size(v::AbstractVariate) = size(v.value)

Base.stride(v::ArrayVariate, k::Int) = stride(v.value, k)


#################### Indexing ####################

Base.getindex(v::ScalarVariate, i::Int) = v.value[i]
Base.getindex(v::ScalarVariate, inds::Union{Range,Vector}) =
  Float64[v[i] for i in inds]
Base.getindex(v::ScalarVariate, ::Colon) = v[[1]]

Base.getindex(v::ArrayVariate, inds::Int...) = getindex(v.value, inds...)

function Base.setindex!(v::ScalarVariate, x, inds)
  length(x) == 1 || throw(ArgumentError("value to set is not of length 1"))
  for i in inds
    i == 1 || throw(BoundsError())
    v.value = x[i]
  end
end
Base.setindex!(v::ScalarVariate, x, ::Colon) = (v[1] = x)

Base.setindex!(v::ArrayVariate, x, inds::Int...) = setindex!(v.value, x, inds...)


#################### I/O ####################

function Base.show(io::IO, v::AbstractVariate)
  print(io, "Object of type \"$(summary(v))\"\n")
  show(io, v.value)
end

function Base.showcompact(io::IO, v::AbstractVariate)
  showcompact(io, v.value)
end


#################### Auxiliary Functions ####################

function names(v::ScalarVariate, prefix)
  AbstractString[string(prefix)]
end

function names(v::ArrayVariate, prefix)
  offset = ndims(v) > 1 ? 1 : 2
  values = similar(v.value, AbstractString)
  for i in 1:length(v)
    s = string(ind2sub(size(v), i))
    values[i] = string(prefix, "[", s[2:end-offset], "]")
  end
  values
end


#################### Mathematical Operators ####################

const BinaryScalarMethods = [
  :(Base.(:(+))),
  :(Base.(:(-))),
  :(Base.(:(*))),
  :(Base.(:(/))),
  :(Base.(:(\))),
  :(Base.(:(^))),
  :(Base.(:(==))),
  :(Base.(:(!=))),
  :(Base.(:(<))),
  :(Base.(:(<=))),
  :(Base.(:(>))),
  :(Base.(:(>=))),
  :(Base.cld),
  :(Base.div),
  :(Base.divrem),
  :(Base.fld),
  :(Base.mod),
  :(Base.rem)
]

for op in BinaryScalarMethods
  @eval ($op)(x::ScalarVariate, y::ScalarVariate) = ($op)(x.value, y.value)
end

const RoundScalarMethods = [
  :(Base.ceil),
  :(Base.floor),
  :(Base.round),
  :(Base.trunc)
]

for op in RoundScalarMethods
  @eval ($op)(x::ScalarVariate) = ($op)(x.value)
  @eval ($op){T}(::Type{T}, x::ScalarVariate) = ($op)(T, x.value)
end

const UnaryScalarMethods = [
  :(Base.(:+)),
  :(Base.(:-)),
  :(Base.abs),
  :(Base.isfinite),
  :(Base.isinf),
  :(Base.isinteger),
  :(Base.isnan),
  :(Base.mod2pi),
  :(Base.one),
  :(Base.sign),
  :(Base.zero)
]

for op in UnaryScalarMethods
  @eval ($op)(x::ScalarVariate) = ($op)(x.value)
end
