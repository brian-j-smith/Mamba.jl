######################################################################
# Distributions Package UnivariateDistribution
######################################################################

#################### Categorical ####################

Distributions.Categorical(p::VectorVariate) =
  Categorical(convert(Vector{Float64}, p))


######################################################################
# Distributions Package MultivariateDistribution
######################################################################

#################### Dirichlet ####################

Dirichlet(alpha::VectorVariate) =
  Dirichlet(convert(Vector{Float64}, alpha))


#################### Multinomial ####################

Multinomial(n::Integer, p::VectorVariate) =
  Multinomial(n, convert(Vector{Float64}, p))


#################### MvNormal ####################

MvNormal(μ::VectorVariate, Σ::Matrix{Float64}) =
  MvNormal(convert(Vector{Float64}, μ), Σ)
MvNormal(μ::Vector{Float64}, Σ::MatrixVariate) =
  MvNormal(μ, convert(Matrix{Float64}, Σ))
MvNormal(μ::VectorVariate, Σ::MatrixVariate) =
  MvNormal(convert(Vector{Float64}, μ), convert(Matrix{Float64}, Σ))

MvNormal(μ::VectorVariate, σ::Vector{Float64}) =
  MvNormal(convert(Vector{Float64}, μ), σ)
MvNormal(μ::Vector{Float64}, σ::VectorVariate) =
  MvNormal(μ, convert(Vector{Float64}, σ))
MvNormal(μ::VectorVariate, σ::VectorVariate) =
  MvNormal(convert(Vector{Float64}, μ), convert(Vector{Float64}, σ))

MvNormal(μ::VectorVariate, σ::Real) =
  MvNormal(convert(Vector{Float64}, μ), σ)

MvNormal(Σ::MatrixVariate) =
  MvNormal(convert(Matrix{Float64}, Σ))

MvNormal(σ::VectorVariate) =
  MvNormal(convert(Vector{Float64}, σ))


#################### MvNormalCanon ####################

MvNormalCanon(h::VectorVariate, J::Matrix{Float64}) =
  MvNormalCanon(convert(Vector{Float64}, h), J)
MvNormalCanon(h::Vector{Float64}, J::MatrixVariate) =
  MvNormalCanon(h, convert(Matrix{Float64}, J))
MvNormalCanon(h::VectorVariate, J::MatrixVariate) =
  MvNormalCanon(convert(Vector{Float64}, h), convert(Matrix{Float64}, J))

MvNormalCanon(h::VectorVariate, prec::Vector{Float64}) =
  MvNormalCanon(convert(Vector{Float64}, h), prec)
MvNormalCanon(h::Vector{Float64}, prec::VectorVariate) =
  MvNormalCanon(h, convert(Vector{Float64}, prec))
MvNormalCanon(h::VectorVariate, prec::VectorVariate) =
  MvNormalCanon(convert(Vector{Float64}, h), convert(Vector{Float64}, prec))

MvNormalCanon(h::VectorVariate, prec::Float64) =
  MvNormalCanon(convert(Vector{Float64}, h), prec)
MvNormalCanon(h::Vector{Float64}, prec::ScalarVariate) =
  MvNormalCanon(h, convert(Float64, prec))
MvNormalCanon(h::VectorVariate, prec::ScalarVariate) =
  MvNormalCanon(convert(Vector{Float64}, h), convert(Float64, prec))

MvNormalCanon(J::MatrixVariate) =
  MvNormalCanon(convert(Matrix{Float64}, J))

MvNormalCanon(prec::VectorVariate) =
  MvNormalCanon(convert(Vector{Float64}, prec))

MvNormalCanon(d::Int, prec::ScalarVariate) =
  MvNormalCanon(d, convert(Float64, prec))


#################### MvTDist ####################

MvTDist(df::ScalarVariate, μ::Vector{Float64}, C::PDMat) =
  MvTDist(convert(Float64, df), μ, C)
MvTDist(df::Float64, μ::VectorVariate, C::PDMat) =
  MvTDist(df, convert(Vector{Float64}, μ), C)
MvTDist(df::ScalarVariate, μ::VectorVariate, C::PDMat) =
  MvTDist(convert(Float64, df), convert(Vector{Float64}, C))

MvTDist(df::ScalarVariate, C::PDMat) =
  MvTDist(convert(Float64, df), C)

MvTDist(df::ScalarVariate, μ::Vector{Float64}, Σ::Matrix{Float64}) =
  MvTDist(convert(Float64, df), μ, Σ)
MvTDist(df::Float64, μ::VectorVariate, Σ::Matrix{Float64}) =
  MvTDist(df, convert(Vector{Float64}, μ), Σ)
MvTDist(df::Float64, μ::Vector{Float64}, Σ::MatrixVariate) =
  MvTDist(df, μ, convert(Matrix{Float64}, Σ))
MvTDist(df::ScalarVariate, μ::VectorVariate, Σ::Matrix{Float64}) =
  MvTDist(convert(Float64, df), convert(Vector{Float64}, μ), Σ)
MvTDist(df::ScalarVariate, μ::Vector{Float64}, Σ::MatrixVariate) =
  MvTDist(convert(Float64, df), μ, convert(Matrix{Float64}, Σ))
MvTDist(df::Float64, μ::VectorVariate, Σ::MatrixVariate) =
  MvTDist(df, convert(Vector{Float64}, μ), convert(Matrix{Float64}, Σ))
MvTDist(df::ScalarVariate, μ::VectorVariate, Σ::MatrixVariate) =
  MvTDist(convert(Float64, df), convert(Vector{Float64}, μ), convert(Matrix{Float64}, Σ))

MvTDist(df::ScalarVariate, Σ::Matrix{Float64}) =
  MvTDist(convert(Float64, df), Σ)
MvTDist(df::Float64, Σ::MatrixVariate) =
  MvTDist(df, convert(Matrix{Float64}, Σ))
MvTDist(df::ScalarVariate, Σ::MatrixVariate) =
  MvTDist(convert(Float64, df), convert(Matrix{Float64}, Σ))


#################### VonMisesFisher ####################

VonMisesFisher(μ::VectorVariate, κ::Real) =
  VonMisesFisher(convert(Vector{Float64}, μ), κ)


######################################################################
# Distributions Package MatrixDistribution
######################################################################

#################### InverseWishart ####################

InverseWishart(df::ScalarVariate, Ψ::Matrix{Float64}) =
  InverseWishart(convert(Float64, df), Ψ)
InverseWishart(df::Real, Ψ::MatrixVariate) =
  InverseWishart(df, convert(Matrix{Float64}, Ψ))
InverseWishart(df::ScalarVariate, Ψ::MatrixVariate) =
  InverseWishart(convert(Float64, df), convert(Matrix{Float64}, Ψ))


#################### Wishart ####################

Wishart(df::ScalarVariate, S::Matrix{Float64}) =
  Wishart(convert(Float64, df), S)
Wishart(df::Real, S::MatrixVariate) =
  Wishart(df, convert(Matrix{Float64}, S))
Wishart(df::ScalarVariate, S::MatrixVariate) =
  Wishart(convert(Float64, df), convert(Matrix{Float64}, S))
