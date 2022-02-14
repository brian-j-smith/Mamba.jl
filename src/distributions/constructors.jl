######################################################################
# Distributions Package UnivariateDistribution
######################################################################

#################### Categorical ####################

Categorical(p::T; check_args = true) where {T <: ArrayVariate{1}} =
  Categorical(convert(Vector{Float64}, p), check_args = check_args)


######################################################################
# Distributions Package MultivariateDistribution
######################################################################

#################### Multinomial ####################

Multinomial(n::Real, p::AbstractVector{T}) where {T<:Real} =
  Multinomial(convert(Int, n), convert(Vector{Float64}, p))


#################### MvNormal ####################

MvNormal(Σ::T) where {T<:ArrayVariate{2}} =
  MvNormal(convert(Matrix{Float64}, Σ))

MvNormal(σ::T) where {T<:ArrayVariate{1}} =
  MvNormal(convert(Vector{Float64}, σ))

PDMat(mat::ArrayStochastic{2}) = PDMat(convert(Matrix{Float64}, mat))


#################### MvNormalCanon ####################

MvNormalCanon(h::AbstractVector{T}, J::U) where {T<:Real, U<:AbstractMatrix} =
  MvNormalCanon(convert(Vector{Float64}, h), convert(Matrix{Float64}, J))

MvNormalCanon(h::AbstractVector{T}, prec::U) where {T<:Real, U<:AbstractVector} =
  MvNormalCanon(convert(Vector{Float64}, h), convert(Vector{Float64}, prec))

MvNormalCanon(prec::AbstractVector{T}) where {T<:Real} =
  MvNormalCanon(convert(Vector{Float64}, prec))

MvNormalCanon(d::Real, prec::T) where {T<:Real} =
  MvNormalCanon(convert(Int, d), convert(Float64, prec))

#################### MvTDist ####################

MvTDist(df::Real, μ::AbstractVector{T}, C::PDMat) where {T<:Real} =
  MvTDist(convert(Float64, df), convert(Vector{Float64}, μ), C)

MvTDist(df::Real, μ::AbstractVector{T}, Σ::AbstractMatrix{U}) where {T<:Real, U<:Real} = 
  MvTDist(convert(Float64, df), convert(Vector{Float64}, μ),
          convert(Matrix{Float64}, Σ))

MvTDist(df::Real, Σ::AbstractMatrix{T}) where {T<:Real} =
  MvTDist(convert(Float64, df), convert(Matrix{Float64}, Σ))


#################### VonMisesFisher ####################

VonMisesFisher(μ::AbstractVector{T}, κ::Real) where {T<:Real} =
  VonMisesFisher(convert(Vector{Float64}, μ), convert(Float64, κ))


######################################################################
# Distributions Package MatrixDistribution
######################################################################

#################### InverseWishart ####################

InverseWishart(df::Real, Ψ::AbstractMatrix{T}) where {T<:Real} =
  InverseWishart(convert(Float64, df), convert(Matrix{Float64}, Ψ))


#################### Wishart ####################

Wishart(df::Real, S::AbstractMatrix{T}) where {T<:Real} =
  Wishart(convert(Float64, df), convert(Matrix{Float64}, S))
