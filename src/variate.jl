#################### Variate Methods ####################

const ArithOps = [
  :(Base.(:+)),
  :(Base.(:.+)),
  :(Base.(:-)),
  :(Base.(:.-)),
  :(Base.(:*)),
  :(Base.(:.*)),
  :(Base.(:/)),
  :(Base.(:./)),
  :(Base.(:.^)),
  :(Base.div),
  :(Base.fld),
  :(Base.mod),
  :(Base.rem)
]

const CompareOps = [
  :(Base.(:(==))),
  :(Base.(:(.==))),
  :(Base.(:(!=))),
  :(Base.(:(.!=))),
  :(Base.(:(>))),
  :(Base.(:(.>))),
  :(Base.(:(>=))),
  :(Base.(:(.>=))),
  :(Base.(:(<))),
  :(Base.(:(.<))),
  :(Base.(:(<=))),
  :(Base.(:(.<=)))
]

const MathOps = [
  :(Base.abs),
  :(Base.acos),
  :(Base.acosh),
  :(Base.asin),
  :(Base.asinh),
  :(Base.atan),
  :(Base.atanh),
  :(Base.ceil),
  :(Base.conj),
  :(Base.cos),
  :(Base.cosh),
  :(Base.digamma),
  :(Base.erf),
  :(Base.erfc),
  :(Base.exp),
  :(Base.exp2),
  :(Base.expm1),
  :(Base.exponent),
  :(Base.floor),
  :(Base.gamma),
  :(Base.lgamma),
  :(Base.log),
  :(Base.log10),
  :(Base.log1p),
  :(Base.log2),
  :(Base.round),
  :(Base.sign),
  :(Base.sin),
  :(Base.sinh),
  :(Base.sqrt),
  :(Base.tan),
  :(Base.tanh),
  :(Base.trunc)
]

const Math2Ops = [
  :(Base.ceil),
  :(Base.floor),
  :(Base.round),
  :(Base.trunc)
]

const UnaryOps = [
  :(Base.(:+)),
  :(Base.(:-)),
  :(Base.ctranspose),
  :(Base.endof),
  :(Base.length),
  :(Base.ndims),
  :(Base.size),
  :(Base.transpose)
]

const VectorOps = [
  :(Base.cummax),
  :(Base.cummin),
  :(Base.cumprod),
  :(Base.cumsum),
  :(Base.cumsum_kbn),
  :(Base.diff),
  :(Base.maximum),
  :(Base.mean),
  :(Base.median),
  :(Base.minimum),
  :(Base.norm),
  :(Base.prod),
  :(Base.std),
  :(Base.sum),
  :(Base.var),
  :(StatsBase.kurtosis),
  :(StatsBase.mad),
  :(StatsBase.skewness)
]

const Vector2Ops = [
  :(Base.cor),
  :(Base.cov),
  :(Base.dot),
  :(StatsBase.corspearman)
]

for op in [ArithOps, CompareOps, Vector2Ops]
  @eval begin
    ($op)(A::Union(Array, Number, SparseMatrixCSC), v::Variate) = ($op)(A, v.data)
    ($op)(v::Variate, A::Union(Array, Number, SparseMatrixCSC)) = ($op)(v.data, A)
    ($op)(v1::Variate, v2::Variate) = ($op)(v1.data, v2.data)
  end
end

for op in [MathOps, UnaryOps, VectorOps]
  @eval ($op)(v::Variate) = ($op)(v.data)
end

for op in Math2Ops
  @eval ($op)(v::Variate, digits::Integer) = ($op)(v.data, digits)
end


function Base.convert{T,N}(::Type{Array{T,N}}, v::Multivariate)
  convert(Array{T,N}, v.data)
end

function Base.convert{T<:Number}(::Type{T}, v::Univariate)
  convert(T, v.data)
end

function Base.getindex(v::Multivariate, inds...)
  getindex(v.data, inds...)
end

function Base.getindex(v::Univariate, inds)
  isa(inds, Number) ? v.data[inds] : map(i -> v.data[i], inds)
end

function Base.setindex!(v::Multivariate, x, inds...)
  setindex!(v.data, x, inds...)
end

function Base.setindex!(v::Univariate, x, inds)
  length(x) == 1 || throw(ErrorException("argument dimensions must match"))
  (length(inds) == 1 && collect(inds)[1] == 1) || throw(BoundsError())
  v.data = x[1]
end

function Base.show(io::IO, v::Variate)
  print(io, "Object of type \"$(summary(v))\"\n")
  show(io, v.data)
  print(io, "\n")
end

function Base.showcompact(io::IO, v::Variate)
  showcompact(io, v.data)
end
