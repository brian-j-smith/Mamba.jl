#################### Flat Distribution ####################

type Flat <: ContinuousUnivariateDistribution end

minimum(d::Flat) = -Inf
maximum(d::Flat) = Inf
insupport(d::Flat, x::Real) = true
insupport{T<:Real}(d::Flat, x::Array{T}) = true

logpdf(d::Flat, x::Real) = 0.0
logpdf{T<:Real}(d::Flat, x::Array{T}) = 0.0

function Truncated(d::Flat, l::Real, u::Real)
  isinf(l) || isinf(u) || error("Flat distributions cannot be bounded")
  Truncated{Flat,Continuous}(d, l, u, 1.0)
end

function insupport{T<:Real}(d::Truncated{Flat}, x::Array{T})
  all(d.lower .<= x .<= d.upper) && insupport(d.untruncated, x)
end
