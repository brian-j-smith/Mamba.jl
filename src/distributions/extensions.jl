#################### Flat Distribution ####################

immutable Flat <: ContinuousUnivariateDistribution end

minimum(d::Flat) = -Inf
maximum(d::Flat) = Inf
insupport(d::Flat, x::Real) = true

logpdf(d::Flat, x::Real) = 0.0

function Truncated(d::Flat, l::Float64, u::Float64)
  Truncated{Flat, Continuous}(d, l, u, 0.0, 1.0, 1.0, 0.0)
end


#################### BDiagNormal Distribution ####################

typealias BDiagNormal MvNormal{PBDiagMat, Vector{Float64}}

function BDiagNormal(μ::Vector{Float64}, Σ::Matrix{Float64})
  n = div(length(μ), size(Σ, 1))
  MvNormal(μ, PBDiagMat(Σ, n))
end

function BDiagNormal(μ::Vector{Float64}, Σ::Vector{Matrix{Float64}})
  MvNormal(μ, PBDiagMat(Σ))
end

BDiagNormal(μ::AbstractVector, Σ::AbstractMatrix) =
  BDiagNormal(convert(Vector{Float64}, μ), convert(Matrix{Float64}, Σ))

BDiagNormal(μ::AbstractVector, Σ::AbstractVector) =
  BDiagNormal(convert(Vector{Float64}, μ), Matrix{Float64}[Σ...])


#################### Null Distribution ####################

immutable NullUnivariateDistribution <: UnivariateDistribution{ValueSupport} end
