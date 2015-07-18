######################################################################
# Univariate Distributions
######################################################################

#################### Categorical ####################

Distributions.Categorical(p::VectorVariate) =
  Categorical(convert(Vector{Float64}, p))


######################################################################
# Multivariate Distributions
######################################################################

#################### BDiagNormal ####################

BDiagNormal(μ::VectorVariate, Σ::Matrix{Float64}) =
  BDiagNormal(convert(Vector{Float64}, μ), Σ)
BDiagNormal(μ::Vector{Float64}, Σ::MatrixVariate) =
  BDiagNormal(μ, convert(Matrix{Float64}, Σ))
BDiagNormal(μ::VectorVariate, Σ::MatrixVariate) =
  BDiagNormal(convert(Vector{Float64}, μ), convert(Matrix{Float64}, Σ))

BDiagNormal(μ, Σ::Vector) =
  BDiagNormal(convert(Vector{Float64}, μ), Matrix{Float64}[Σ...])


#################### Dirichlet ####################

Distributions.Dirichlet(alpha::VectorVariate) =
  Dirichlet(convert(Vector{Float64}, alpha))


#################### Multinomial ####################

Distributions.Multinomial(n::Integer, p::VectorVariate) =
  Multinomial(n, convert(Vector{Float64}, p))


#################### MvNormal ####################

Distributions.MvNormal(μ::VectorVariate, Σ::Matrix{Float64}) =
  MvNormal(convert(Vector{Float64}, μ), Σ)
Distributions.MvNormal(μ::Vector{Float64}, Σ::MatrixVariate) =
  MvNormal(μ, convert(Matrix{Float64}, Σ))
Distributions.MvNormal(μ::VectorVariate, Σ::MatrixVariate) =
  MvNormal(convert(Vector{Float64}, μ), convert(Matrix{Float64}, Σ))

Distributions.MvNormal(μ::VectorVariate, σ::Vector{Float64}) =
  MvNormal(convert(Vector{Float64}, μ), σ)
Distributions.MvNormal(μ::Vector{Float64}, σ::VectorVariate) =
  MvNormal(μ, convert(Vector{Float64}, σ))
Distributions.MvNormal(μ::VectorVariate, σ::VectorVariate) =
  MvNormal(convert(Vector{Float64}, μ), convert(Vector{Float64}, σ))

Distributions.MvNormal(μ::VectorVariate, σ::Real) =
  MvNormal(convert(Vector{Float64}, μ), σ)

Distributions.MvNormal(Σ::MatrixVariate) =
  MvNormal(convert(Matrix{Float64}, Σ))

Distributions.MvNormal(σ::VectorVariate) =
  MvNormal(convert(Vector{Float64}, σ))


#################### MvNormalCanon ####################

Distributions.MvNormalCanon(h::VectorVariate, J::Matrix{Float64}) =
  MvNormalCanon(convert(Vector{Float64}, h), J)
Distributions.MvNormalCanon(h::Vector{Float64}, J::MatrixVariate) =
  MvNormalCanon(h, convert(Matrix{Float64}, J))
Distributions.MvNormalCanon(h::VectorVariate, J::MatrixVariate) =
  MvNormalCanon(convert(Vector{Float64}, h), convert(Matrix{Float64}, J))

Distributions.MvNormalCanon(h::VectorVariate, prec::Vector{Float64}) =
  MvNormalCanon(convert(Vector{Float64}, h), prec)
Distributions.MvNormalCanon(h::Vector{Float64}, prec::VectorVariate) =
  MvNormalCanon(h, convert(Vector{Float64}, prec))
Distributions.MvNormalCanon(h::VectorVariate, prec::VectorVariate) =
  MvNormalCanon(convert(Vector{Float64}, h), convert(Vector{Float64}, prec))

Distributions.MvNormalCanon(h::VectorVariate, prec::Float64) =
  MvNormalCanon(convert(Vector{Float64}, h), prec)
Distributions.MvNormalCanon(h::Vector{Float64}, prec::ScalarVariate) =
  MvNormalCanon(h, convert(Float64, prec))
Distributions.MvNormalCanon(h::VectorVariate, prec::ScalarVariate) =
  MvNormalCanon(convert(Vector{Float64}, h), convert(Float64, prec))

Distributions.MvNormalCanon(J::MatrixVariate) =
  MvNormalCanon(convert(Matrix{Float64}, J))

Distributions.MvNormalCanon(prec::VectorVariate) =
  MvNormalCanon(convert(Vector{Float64}, prec))

Distributions.MvNormalCanon(d::Int, prec::ScalarVariate) =
  MvNormalCanon(d, convert(Float64, prec))


#################### MvTDist ####################

Distributions.MvTDist(df::ScalarVariate, μ::Vector{Float64}, C::PDMat) =
  MvTDist(convert(Float64, df), μ, C)
Distributions.MvTDist(df::Float64, μ::VectorVariate, C::PDMat) =
  MvTDist(df, convert(Vector{Float64}, μ), C)
Distributions.MvTDist(df::ScalarVariate, μ::VectorVariate, C::PDMat) =
  MvTDist(convert(Float64, df), convert(Vector{Float64}, C))

Distributions.MvTDist(df::ScalarVariate, C::PDMat) =
  MvTDist(convert(Float64, df), C)

Distributions.MvTDist(df::ScalarVariate, μ::Vector{Float64}, Σ::Matrix{Float64}) =
  MvTDist(convert(Float64, df), μ, Σ)
Distributions.MvTDist(df::Float64, μ::VectorVariate, Σ::Matrix{Float64}) =
  MvTDist(df, convert(Vector{Float64}, μ), Σ)
Distributions.MvTDist(df::Float64, μ::Vector{Float64}, Σ::MatrixVariate) =
  MvTDist(df, μ, convert(Matrix{Float64}, Σ))
Distributions.MvTDist(df::ScalarVariate, μ::VectorVariate, Σ::Matrix{Float64}) =
  MvTDist(convert(Float64, df), convert(Vector{Float64}, μ), Σ)
Distributions.MvTDist(df::ScalarVariate, μ::Vector{Float64}, Σ::MatrixVariate) =
  MvTDist(convert(Float64, df), μ, convert(Matrix{Float64}, Σ))
Distributions.MvTDist(df::Float64, μ::VectorVariate, Σ::MatrixVariate) =
  MvTDist(df, convert(Vector{Float64}, μ), convert(Matrix{Float64}, Σ))
Distributions.MvTDist(df::ScalarVariate, μ::VectorVariate, Σ::MatrixVariate) =
  MvTDist(convert(Float64, df), convert(Vector{Float64}, μ), convert(Matrix{Float64}, Σ))

Distributions.MvTDist(df::ScalarVariate, Σ::Matrix{Float64}) =
  MvTDist(convert(Float64, df), Σ)
Distributions.MvTDist(df::Float64, Σ::MatrixVariate) =
  MvTDist(df, convert(Matrix{Float64}, Σ))
Distributions.MvTDist(df::ScalarVariate, Σ::MatrixVariate) =
  MvTDist(convert(Float64, df), convert(Matrix{Float64}, Σ))


#################### VonMisesFisher ####################

Distributions.VonMisesFisher(μ::VectorVariate, κ::Real) =
  VonMisesFisher(convert(Vector{Float64}, μ), κ)


######################################################################
# Matrix Distributions
######################################################################

#################### InverseWishart ####################

Distributions.InverseWishart(df::ScalarVariate, Ψ::Matrix{Float64}) =
  InverseWishart(convert(Float64, df), Ψ)
Distributions.InverseWishart(df::Real, Ψ::MatrixVariate) =
  InverseWishart(df, convert(Matrix{Float64}, Ψ))
Distributions.InverseWishart(df::ScalarVariate, Ψ::MatrixVariate) =
  InverseWishart(convert(Float64, df), convert(Matrix{Float64}, Ψ))


#################### Wishart ####################

Distributions.Wishart(df::ScalarVariate, S::Matrix{Float64}) =
  Wishart(convert(Float64, df), S)
Distributions.Wishart(df::Real, S::MatrixVariate) =
  Wishart(df, convert(Matrix{Float64}, S))
Distributions.Wishart(df::ScalarVariate, S::MatrixVariate) =
  Wishart(convert(Float64, df), convert(Matrix{Float64}, S))
