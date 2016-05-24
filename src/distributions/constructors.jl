######################################################################
# Distributions Package UnivariateDistribution
######################################################################

#################### Categorical ####################

Categorical{T<:Real}(p::AbstractVector{T}) =
  Categorical(convert(Vector{Float64}, p))


######################################################################
# Distributions Package MultivariateDistribution
######################################################################

#################### Dirichlet ####################

Dirichlet{T<:Real}(alpha::AbstractVector{T}) =
  Dirichlet(convert(Vector{Float64}, alpha))


#################### Multinomial ####################

Multinomial{T<:Real}(n::Real, p::AbstractVector{T}) =
  Multinomial(convert(Int, n), convert(Vector{Float64}, p))


#################### MvNormal ####################

MvNormal{T<:Real, U<:Real}(μ::AbstractVector{T}, Σ::AbstractMatrix{U}) =
  MvNormal(convert(Vector{Float64}, μ), convert(Matrix{Float64}, Σ))

MvNormal{T<:Real, U<:Real}(μ::AbstractVector{T}, σ::AbstractVector{U}) =
  MvNormal(convert(Vector{Float64}, μ), convert(Vector{Float64}, σ))

MvNormal{T<:Real}(μ::AbstractVector{T}, σ::Real) =
  MvNormal(convert(Vector{Float64}, μ), convert(Float64, σ))

MvNormal{T<:Real}(Σ::AbstractMatrix{T}) =
  MvNormal(convert(Matrix{Float64}, Σ))

MvNormal{T<:Real}(σ::AbstractVector{T}) =
  MvNormal(convert(Vector{Float64}, σ))


#################### MvNormalCanon ####################

MvNormalCanon{T<:Real, U<:AbstractMatrix}(h::AbstractVector{T}, J::U) =
  MvNormalCanon(convert(Vector{Float64}, h), convert(Matrix{Float64}, J))

MvNormalCanon{T<:Real, U<:AbstractVector}(h::AbstractVector{T}, prec::U) =
  MvNormalCanon(convert(Vector{Float64}, h), convert(Vector{Float64}, prec))

MvNormalCanon{T<:Real, U<:Real}(h::AbstractVector{T}, prec::U) =
  MvNormalCanon(convert(Vector{Float64}, h), convert(Float64, prec))

MvNormalCanon{T<:Real}(J::AbstractMatrix{T}) =
  MvNormalCanon(convert(Matrix{Float64}, J))

MvNormalCanon{T<:Real}(prec::AbstractVector{T}) =
  MvNormalCanon(convert(Vector{Float64}, prec))

MvNormalCanon{T<:Real}(d::Real, prec::T) =
  MvNormalCanon(convert(Int, d), convert(Float64, prec))

#################### MvTDist ####################

MvTDist{T<:Real}(df::Real, μ::AbstractVector{T}, C::PDMat) =
  MvTDist(convert(Float64, df), convert(Vector{Float64}, μ), C)

MvTDist{T<:Real, U<:Real}(df::Real, μ::AbstractVector{T}, Σ::AbstractMatrix{U}) =
  MvTDist(convert(Float64, df), convert(Vector{Float64}, μ),
          convert(Matrix{Float64}, Σ))

MvTDist{T<:Real}(df::Real, Σ::AbstractMatrix{T}) =
  MvTDist(convert(Float64, df), convert(Matrix{Float64}, Σ))


#################### VonMisesFisher ####################

VonMisesFisher{T<:Real}(μ::AbstractVector{T}, κ::Real) =
  VonMisesFisher(convert(Vector{Float64}, μ), convert(Float64, κ))


######################################################################
# Distributions Package MatrixDistribution
######################################################################

#################### InverseWishart ####################

InverseWishart{T<:Real}(df::Real, Ψ::AbstractMatrix{T}) =
  InverseWishart(convert(Float64, df), convert(Matrix{Float64}, Ψ))


#################### Wishart ####################

Wishart{T<:Real}(df::Real, S::AbstractMatrix{T}) =
  Wishart(convert(Float64, df), convert(Matrix{Float64}, S))
