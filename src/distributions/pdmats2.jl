#################### PBDiagMat Block-Diagonal Matrices ####################

module PDMats2

  import Base: +, *, /, \
  import Base: diag, full, inv, logdet, size
  import Base.LinAlg: Cholesky
  import PDMats: AbstractPDMat, dim, invquad, invquad!, quad, quad!,
         whiten, whiten!, unwhiten, unwhiten!

  export PBDiagMat


  #################### Types and Constructors ####################

  type PBDiagMat <: AbstractPDMat
    dim::Int
    mat::SparseMatrixCSC{Float64, Int}
    chol::Vector{Cholesky{Float64}}
    scale::Int
  end

  function PBDiagMat{T<:Real}(v::Vector{Matrix{T}}, n::Integer=1)
    mat = spbdiagm(v, n)
    chol = map(cholfact, v)
    PBDiagMat(size(mat, 1), mat, chol, n)
  end

  function PBDiagMat{T<:Real}(x::Matrix{T}, n::Integer=1)
    PBDiagMat(Matrix{T}[x], n)
  end


  #################### Base Methods ####################

  +(a::PBDiagMat, b::Matrix{Float64}) = a.mat + b
  +(a::Matrix{Float64}, b::PBDiagMat) = b + a

  *(a::PBDiagMat, c::Float64) = mapchol(x -> c * full(x), a)
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
  full(a::PBDiagMat) = full(a.mat)
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
      r[j] = sumabs2(wx[:, j])
    end
    r
  end

  invquad(a::PBDiagMat, x::Vector{Float64}) = sumabs2(whiten(a, x))


  #################### Auxiliary Functions ####################

  mapchol(f::Function, a::PBDiagMat) = PBDiagMat(map(f, a.chol), a.scale)

  function spbdiagm{T<:Real}(v::Union{Vector{Matrix{T}},
                                      Vector{UpperTriangular{T, Matrix{T}}}},
                             n::Integer=1)
    vn = [fill(v, n)...;]

    len = mapreduce(splength, +, vn)
    I = Array{Int}(len)
    J = Array{Int}(len)
    V = Array{Float64}(len)

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

  function splength{T}(x::UpperTriangular{T, Matrix{T}})
    m, n = minmax(size(x)...)
    Int(m * (m + 1) / 2) + (n - m) * m
  end
  isnonzero{T}(x::UpperTriangular{T, Matrix{T}}, i::Integer, j::Integer) =
    j >= i

end
