#################### BDiagNormal Distribution ####################

typealias BDiagNormal MvNormal{PBDiagMat,Vector{Float64}}

function BDiagNormal(μ::Vector{Float64}, Σ::Matrix{Float64})
  n = div(length(μ), size(Σ, 1))
  MvNormal(μ, PBDiagMat(Σ, n))
end

function BDiagNormal(μ::Vector{Float64}, Σ::Vector{Matrix{Float64}})
  MvNormal(μ, PBDiagMat(Σ))
end

BDiagNormal(μ::VectorVariate, Σ::Matrix{Float64}) =
  BDiagNormal(convert(Vector{Float64}, μ), Σ)
BDiagNormal(μ::Vector{Float64}, Σ::MatrixVariate) =
  BDiagNormal(μ, convert(Matrix{Float64}, Σ))
BDiagNormal(μ::VectorVariate, Σ::MatrixVariate) =
  BDiagNormal(convert(Vector{Float64}, μ), convert(Matrix{Float64}, Σ))

BDiagNormal(μ, Σ::Vector) =
  BDiagNormal(convert(Vector{Float64}, μ), Matrix{Float64}[Σ...])
