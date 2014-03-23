#################### Null Distribution ####################

type NullDistribution <: Distribution end


#################### Flat Distribution ####################

type Flat <: Distribution end

insupport(d::Flat, x::Real) = true
insupport{T<:Real}(d::Flat, x::VecOrMat{T}) = true
logpdf(d::Flat, x::Real) = 0
logpdf{T<:Real}(d::Flat, x::VecOrMat{T}) = 0
