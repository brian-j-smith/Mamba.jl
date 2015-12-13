## Deprecated at 0.7.2

using Distributions

function Logical(expr::Expr, monitor::Union{Bool, Vector{Int}}=true)
  msg = string("Logical(expr::Expr, monitor::$(typeof(monitor))) is deprecated; ",
               "use Logical(f::Function, monitor::$(typeof(monitor))) instead")
  Base.depwarn(msg, :Logical)

  value = Float64(NaN)
  l = ScalarLogical(value, :nothing, Int[], modelfx(depfxargs, expr),
                    depsrc(expr), Symbol[])
  setmonitor!(l, monitor)
end

function Logical(d::Integer, expr::Expr, monitor::Union{Bool, Vector{Int}}=true)
  msg = string("Logical(d::Integer, expr::Expr, monitor::$(typeof(monitor))) is deprecated; ",
               "use Logical(d::Integer, f::Function, monitor::$(typeof(monitor))) instead")
  Base.depwarn(msg, :Logical)

  value = Array(Float64, fill(0, d)...)
  l = ArrayLogical(value, :nothing, Int[], modelfx(depfxargs, expr),
                   depsrc(expr), Symbol[])
  setmonitor!(l, monitor)
end

function Stochastic(expr::Expr, monitor::Union{Bool, Vector{Int}}=true)
  msg = string("Stochastic(expr::Expr, monitor::$(typeof(monitor))) is deprecated; ",
               "use Stochastic(f::Function, monitor::$(typeof(monitor))) instead")
  Base.depwarn(msg, :Stochastic)

  value = Float64(NaN)
  s = ScalarStochastic(value, :nothing, Int[], modelfx(depfxargs, expr),
                       depsrc(expr), Symbol[], NullUnivariateDistribution())
  setmonitor!(s, monitor)
end

function Stochastic(d::Integer, expr::Expr,
                    monitor::Union{Bool, Vector{Int}}=true)
  msg = string("Stochastic(d::Integer, expr::Expr, monitor::$(typeof(monitor))) is deprecated; ",
               "use Stochastic(d::Integer, f::Function, monitor::$(typeof(monitor))) instead")
  Base.depwarn(msg, :Stochastic)

  value = Array(Float64, fill(0, d)...)
  s = ArrayStochastic(value, :nothing, Int[], modelfx(depfxargs, expr),
                      depsrc(expr), Symbol[], NullUnivariateDistribution())
  setmonitor!(s, monitor)
end

function Sampler(params::Vector{Symbol}, expr::Expr, tune::Dict=Dict())
  msg = string("Sampler(params::Vector{Symbol}, expr::Expr, tune::Dict) is deprecated; ",
               "use Sampler(params::Vector{Symbol}, f::Function, tune::Dict) instead")
  Base.depwarn(msg, :Sampler)

  Sampler(params, modelfx(samplerfxargs, expr), tune, Symbol[])
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

function modelfx(literalargs::Vector{Tuple{Symbol, Symbol}}, expr::Expr)
  args = Expr(:tuple, map(arg -> Expr(:(::), arg[1], arg[2]), literalargs)...)
  eval(Expr(:function, args, expr))
end

function depsrc(expr::Expr)
  if expr.head == :ref && expr.args[1] == :model && isa(expr.args[2], QuoteNode)
    Symbol[expr.args[2].value]
  else
    mapreduce(depsrc, union, expr.args)
  end
end

function depsrc(expr)
  Symbol[]
end


## Deprecated at 0.7.3

export bmmg!, BMMG, BMMGVariate

@deprecate(bmmg!, bmc3!)
@deprecate(BMMG, BMC3)
@deprecate(BMMGVariate, BMC3Variate)

function BMC3(params::Vector{Symbol}, d::Integer, k::Integer=1)
  msg = string("BMC3(params::Vector{Symbol}, d::Integer, k::Integer=1) is deprecated; ",
               "use BMC3(params::Vector{Symbol}; k::Integer=1) instead")
  Base.depwarn(msg, :BMC3)

  BMC3(params, k=k)
end


## Deprecated at 0.7.4

macro depsamplermethod(T)
  esc(quote
        function $T{T<:Real}(x::AbstractVector{T}, tune::Void)
          msg = string($T, "{T<:Real}(x::AbstractVector{T}, tune::Void) is deprecated; ",
                       "use ", $T, "{T<:Real}(x::AbstractVector{T}) instead")
          Base.depwarn(msg, symbol($T))
          $T(x)
        end
      end)
end

@depsamplermethod AMMVariate
@depsamplermethod AMWGVariate
@depsamplermethod BHMCVariate
@depsamplermethod BMC3Variate
@depsamplermethod BMGVariate
@depsamplermethod DGSVariate
@depsamplermethod HMCVariate
@depsamplermethod MALAVariate
@depsamplermethod NUTSVariate
@depsamplermethod SliceVariate
@depsamplermethod SliceSimplexVariate
