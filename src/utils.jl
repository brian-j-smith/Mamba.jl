#################### Model Expression Operators ####################

function modelfx(literalargs::Vector{Tuple{Symbol, Symbol}}, f::Function)
  modelfxsrc(literalargs, f)[1]
end

function modelfxsrc(literalargs::Vector{Tuple{Symbol, Symbol}}, f::Function)
  args = Expr(:tuple, map(arg -> Expr(:(::), arg[1], arg[2]), literalargs)...)
  expr, src = modelexprsrc(f, literalargs)
  fx = eval(Expr(:function, args, expr))
  (fx, src)
end


function modelexprsrc(f::Function, literalargs::Vector{Tuple{Symbol, Symbol}})
  ast = ccall(:jl_uncompress_ast, Any, (Any, Any), f.code, f.code.ast)
  fargs = ast.args[1]
  n = length(fargs)
  fkeys = Symbol[arg.args[1] for arg in fargs]
  ftypes = Symbol[arg.args[2] for arg in fargs]

  literalinds = Int[]
  for (key, T) in literalargs
    i = findfirst(fkey -> fkey == key, fkeys)
    if i != 0 && ftypes[i] == T
      push!(literalinds, i)
    end
  end
  nodeinds = setdiff(1:n, literalinds)

  all(T -> T == :Any, ftypes[nodeinds]) ||
    throw(ArgumentError("model node arguments are not all of type Any"))

  modelargs = Array(Any, n)
  for i in nodeinds
    modelargs[i] = Expr(:ref, :model, QuoteNode(fkeys[i]))
  end
  for i in literalinds
    modelargs[i] = fkeys[i]
  end
  expr = Expr(:block, Expr(:(=), :f, f), Expr(:call, :f, modelargs...))

  (expr, fkeys[nodeinds])
end


#################### Mathematical Operators ####################

isprobvec(p::AbstractVector) = isprobvec(convert(Vector{Float64}, p))

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


#################### Auxiliary Functions ####################

## pmap2 is a partial work-around for the pmap issue in julia 0.4.0 of worker
## node errors being blocked.  In single-processor mode, pmap2 calls map
## instead to avoid the error handling issue.  In multi-processor model, pmap is
## called and will apply its error processing.  If and when the pmap issue is
## resolved in a future version of julia, calls to pmap2 should be reverted to
## a call to pmap.

function pmap2(f::Function, lsts::AbstractArray)
  nprocs() > 1 ? pmap(f, lsts) : map(f, lsts)
end
