#################### Null Distribution ####################

type NullDistribution <: Distribution end


#################### Flat Distribution ####################

immutable Flat <: Distribution
  length::Integer
  Flat(length::Real) = new(integer(length))
end

insupport{T<:Real}(d::Flat, x::Union(T, Vector{T})) = d.length == length(x)
function logpdf{T<:Real}(d::Flat, x::Union(T, Vector{T}))
  d.length == length(x) ? 0 : throw(BoundsError())
end
