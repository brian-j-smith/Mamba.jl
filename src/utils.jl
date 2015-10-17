#################### Model Expression Operators ####################

function modelexpr(f::Function)
  ast = ccall(:jl_uncompress_ast, Any, (Any, Any), f.code, f.code.ast)

  types = Symbol[args.args[2] for args in ast.args[1]]
  all(T -> T == :Any, types) ||
    throw(ArgumentError("node function arguments must all be of type Any"))

  src = Symbol[args.args[1] for args in ast.args[1]]
  modelargs = [Expr(:ref, :model, QuoteNode(s)) for s in src]
  Expr(:block, Expr(:(=), :f, f), Expr(:call, :f, modelargs...))
end

macro modelexpr(args...)
  quote
    n = length($args)
    ex = Array(Expr, n)
    for i in 1:(n-1)
      x = $args[i]
      ex[i] = Expr(:(=), x, Expr(:ref, :model, QuoteNode(x)))
    end
    ex[n] = $args[n]
    Expr(:block, ex...)
  end
end

macro promote_scalarvariate(V)
  quote
    Base.promote_rule{T<:Real}(::Type{$V}, ::Type{T}) = Float64
  end
end


#################### Mathematical Operators ####################

cummean(x::AbstractArray) = mapslices(cummean, x, 1)

function cummean{T<:Real}(x::AbstractVector{T})
  y = similar(x, Float64)
  xs = 0.0
  for i in 1:length(x)
    xs += x[i]
    y[i] = xs / i
  end
  y
end

dot(x) = dot(x, x)

invlogit(x::Real) = 1.0 / (exp(-x) + 1.0)
invlogit(x::AbstractArray) = map(invlogit, x)

logit(x::Real) = log(x / (1.0 - x))
logit(x::AbstractArray) = map(logit, x)

## Csorgo S and Faraway JJ. The exact and asymptotic distributions of the
## Cramer-von Mises statistic. Journal of the Royal Statistical Society,
## Series B, 58: 221-234, 1996.
function pcramer(q::Real)
  p = 0.0
  for k in 0:3
    c1 = 4.0 * k + 1.0
    c2 = c1^2 / (16.0 * q)
    p += gamma(k + 0.5) / factorial(k) * sqrt(c1) * exp(-c2) * besselk(0.25, c2)
  end
  p / (pi^1.5 * sqrt(q))
end
