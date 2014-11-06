#################### BDiagNormal ####################

typealias BDiagNormal MvNormal{PBDiagMat,Vector{Float64}}

function BDiagNormal(μ::Vector{Float64}, Σ::Matrix{Float64})
  n = div(length(μ), size(Σ, 1))
  MvNormal(μ, PBDiagMat(Σ, n))
end

function BDiagNormal(μ::Vector{Float64}, Σ::Vector{Matrix{Float64}})
  MvNormal(μ, PBDiagMat(Σ))
end
