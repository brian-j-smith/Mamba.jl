######################################################################
# Distributions Package UnivariateDistribution
######################################################################

#################### Categorical ####################

Distributions.Categorical(p::AbstractVector) =
  Categorical(convert(Vector{Float64}, p))


######################################################################
# Distributions Package MultivariateDistribution
######################################################################

#################### Dirichlet ####################

Dirichlet(alpha::AbstractVector) =
  Dirichlet(convert(Vector{Float64}, alpha))


#################### Multinomial ####################

Multinomial(n::Real, p::AbstractVector) =
  Multinomial(convert(Int, n), convert(Vector{Float64}, p))


#################### MvNormal ####################

MvNormal(μ::AbstractVector, Σ::AbstractMatrix) =
  MvNormal(convert(Vector{Float64}, μ), convert(Matrix{Float64}, Σ))

MvNormal(μ::AbstractVector, σ::AbstractVector) =
  MvNormal(convert(Vector{Float64}, μ), convert(Vector{Float64}, σ))

MvNormal(μ::AbstractVector, σ::Real) =
  MvNormal(convert(Vector{Float64}, μ), convert(Float64, σ))

MvNormal(Σ::AbstractMatrix) =
  MvNormal(convert(Matrix{Float64}, Σ))

MvNormal(σ::AbstractVector) =
  MvNormal(convert(Vector{Float64}, σ))


#################### MvNormalCanon ####################

MvNormalCanon(h::AbstractVector, J::AbstractMatrix) =
  MvNormalCanon(convert(Vector{Float64}, h), convert(Matrix{Float64}, J))

MvNormalCanon(h::AbstractVector, prec::AbstractVector) =
  MvNormalCanon(convert(Vector{Float64}, h), convert(Vector{Float64}, prec))

MvNormalCanon(h::AbstractVector, prec::Real) =
  MvNormalCanon(convert(Vector{Float64}, h), convert(Float64, prec))

MvNormalCanon(J::AbstractMatrix) =
  MvNormalCanon(convert(Matrix{Float64}, J))

MvNormalCanon(prec::AbstractVector) =
  MvNormalCanon(convert(Vector{Float64}, prec))

MvNormalCanon(d::Real, prec::Real) =
  MvNormalCanon(convert(Int, d), convert(Float64, prec))


#################### MvTDist ####################

MvTDist(df::Real, μ::AbstractVector, C::PDMat) =
  MvTDist(convert(Float64, df), convert(Vector{Float64}, μ), C)

MvTDist(df::Real, C::PDMat) =
  MvTDist(convert(Float64, df), C)

MvTDist(df::Real, μ::AbstractVector, Σ::AbstractMatrix) =
  MvTDist(convert(Float64, df), convert(Vector{Float64}, μ),
          convert(Matrix{Float64}, Σ))

MvTDist(df::Real, Σ::AbstractMatrix) =
  MvTDist(convert(Float64, df), convert(Matrix{Float64}, Σ))


#################### VonMisesFisher ####################

VonMisesFisher(μ::AbstractVector, κ::Real) =
  VonMisesFisher(convert(Vector{Float64}, μ), convert(Float64, κ))


######################################################################
# Distributions Package MatrixDistribution
######################################################################

#################### InverseWishart ####################

InverseWishart(df::Real, Ψ::AbstractMatrix) =
  InverseWishart(convert(Float64, df), convert(Matrix{Float64}, Ψ))


#################### Wishart ####################

Wishart(df::Real, S::AbstractMatrix) =
  Wishart(convert(Float64, df), convert(Matrix{Float64}, S))
