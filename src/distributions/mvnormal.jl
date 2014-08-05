#################### BDiagNormal ####################

typealias BDiagNormal GenericMvNormal{PBDiagMat}

function BDiagNormal(μ::Vector{Float64}, Σ::Matrix{Float64})
  n = div(length(μ), size(Σ, 1))
  GenericMvNormal(μ, PBDiagMat(Σ, n))
end

function BDiagNormal(μ::Vector{Float64}, Σ::Vector{Matrix{Float64}})
  GenericMvNormal(μ, PBDiagMat(Σ))
end
