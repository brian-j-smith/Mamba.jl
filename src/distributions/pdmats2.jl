#################### PBDiagMat Block-Diagonal Matrices ####################

module PDMats2

  using SparseArrays: SparseMatrixCSC
  using LinearAlgebra: Cholesky, UpperTriangular
  using PDMats: AbstractPDMat

  import Base: +, *, /, \, size
  import LinearAlgebra: diag, inv, logdet
  import PDMats: dim, invquad, invquad!, quad, quad!,
         whiten, whiten!, unwhiten, unwhiten!

  export PBDiagMat

  #################### Types and Constructors ####################

  struct PBDiagMat <: AbstractPDMat{Float64}
    dim::Int
    mat::SparseMatrixCSC{Float64, Int}
    chol::Vector{Cholesky{Float64}}
    scale::Int
  end

  function PBDiagMat(v::Vector{Matrix{T}}, n::Integer=1) where {T<:Real}
    mat = spbdiagm(v, n)
    chol = map(cholesky, v)
    PBDiagMat(size(mat, 1), mat, chol, n)
  end

  function PBDiagMat(x::Matrix{T}, n::Integer=1) where {T<:Real}
    PBDiagMat(Matrix{T}[x], n)
  end

  #################### Base Methods ####################

  +(a::PBDiagMat, b::Matrix{Float64}) = a.mat + b
  +(a::Matrix{Float64}, b::PBDiagMat) = b + a

  *(a::PBDiagMat, c::Float64) = mapchol(x -> c * Matrix(x), a)
  /(a::PBDiagMat, c::Float64) = a * inv(c)
  *(c::Float64, a::PBDiagMat) = a * c

  *(a::PBDiagMat, x::StridedVecOrMat) = a.mat * x
  *(x::StridedVecOrMat, a::PBDiagMat) = x * a.mat

  \(a::PBDiagMat, x::StridedVecOrMat) = inv(a) * x
  \(a::Matrix{Float64}, b::PBDiagMat) = inv(a) * b

  dim(a::PBDiagMat) = a.dim
  size(a::PBDiagMat) = size(a.mat)
  size(a::PBDiagMat, i) = size(a)[i]

  diag(a::PBDiagMat) = diag(a.mat)
  Matrix(a::PBDiagMat) = Matrix(a.mat)
  inv(a::PBDiagMat) = mapchol(inv, a)
  logdet(a::PBDiagMat) = a.scale * mapreduce(logdet, +, a.chol)


  #################### Whiten and Unwhiten ####################

  function whiten(a::PBDiagMat, x::DenseVecOrMat{Float64})
    au_inv = map(ac -> inv(ac[:U]), a.chol)
    At_mul_B(spbdiagm(au_inv, a.scale), x)
  end

  function whiten!(a::PBDiagMat, x::DenseVecOrMat{Float64})
    x[:] = whiten(a, x)
  end

  function unwhiten(a::PBDiagMat, x::DenseVecOrMat{Float64})
    au = map(ac -> ac[:U], a.chol)
    At_mul_B(spbdiagm(au, a.scale), x)
  end

  function unwhiten!(a::PBDiagMat, x::DenseVecOrMat{Float64})
    x[:] = unwhiten(a, x)
  end


  #################### Quadratic Forms ####################

  function quad!(r::Array{Float64}, a::PBDiagMat, x::Matrix{Float64})
    n = size(x, 2)
    length(r) == n || throw(ArgumentError("incompatible array length"))
    ax = a * x
    for j = 1:n
      r[j] = dot(a[:, j], ax[:, j])
    end
    r
  end

  quad(a::PBDiagMat, x::Vector{Float64}) = dot(x, a * x)

  function invquad!(r::Array{Float64}, a::PBDiagMat, x::Matrix{Float64})
    n = size(x, 2)
    length(r) == n || throw(ArgumentError("incompatible array length"))
    wx = whiten(a, x)
    for j = 1:n
      r[j] = sum(abs2, wx[:, j])
    end
    r
  end

  invquad(a::PBDiagMat, x::Vector{Float64}) = sum(abs2, whiten(a, x))


  #################### Auxiliary Functions ####################

  mapchol(f::Function, a::PBDiagMat) = PBDiagMat(map(f, a.chol), a.scale)

  const VecOfMatsOrUppers{T} = Union{Vector{Matrix{T}}, Vector{UpperTriangular{T, Matrix{T}}}}
  function spbdiagm(v::VecOfMatsOrUppers{T}, n::Integer=1) where {T<:Real}
    vn = [fill(v, n)...;]

    len = mapreduce(splength, +, vn)
    I = Array{Int}(undef, len)
    J = Array{Int}(undef, len)
    V = Array{Float64}(undef, len)

    k = 1
    offset = 0
    for x in vn
      m = size(x, 1)
      size(x, 2) == m || throw(ArgumentError("blocks are not square matrices"))
      for i in 1:m, j in 1:m
        if isnonzero(x, i, j)
          I[k] = offset + i
          J[k] = offset + j
          V[k] = x[i, j]
          k += 1
        end
      end
      offset += m
    end

    sparse(I, J, V, offset, offset)
  end

  splength(x::Matrix) = length(x)
  isnonzero(x::Matrix, i::Integer, j::Integer) = true

  function splength(x::UpperTriangular{T, Matrix{T}}) where {T}
    m, n = minmax(size(x)...)
    Int(m * (m + 1) / 2) + (n - m) * m
  end
  isnonzero(x::UpperTriangular{T, Matrix{T}}, i::Integer, j::Integer) where {T} =
    j >= i

end
