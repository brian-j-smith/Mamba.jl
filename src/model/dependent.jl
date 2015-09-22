#################### Dependent Methods ####################

function Base.show(io::IO, d::AbstractDependent)
  msg = string(ifelse(length(d.monitor) > 0, "A ", "An un"),
               "monitored node of type \"", summary(d), "\"\n")
  print(io, msg)
  show(io, d.value)
end

function Base.showall(io::IO, d::AbstractDependent)
  show(io, d)
  print(io, "\nFunction:\n")
  show(io, d.eval.code)
  print(io, "\n\nSource Nodes:\n")
  show(io, d.sources)
  print(io, "\n\nTarget Nodes:\n")
  show(io, d.targets)
end

dims(d::AbstractDependent) = size(d)

invlink(d::AbstractDependent, x, transform::Bool=true) = x

link(d::AbstractDependent, x, transform::Bool=true) = x

logpdf(d::AbstractDependent, transform::Bool=false) = 0.0

function names(d::AbstractDependent)
  names(d, d.symbol)
end

function setmonitor!(d::AbstractDependent, monitor::Bool)
  value = monitor ? Int[0] : Int[]
  setmonitor!(d, value)
end

function setmonitor!(d::AbstractDependent, monitor::Vector{Int})
  values = monitor
  d.linklength = length(link(d, d.value, false))
  if d.linklength > 0 && length(monitor) > 0
    if monitor[1] == 0
      values = collect(1:d.linklength)
    elseif minimum(monitor) < 1 || maximum(monitor) > d.linklength
      throw(BoundsError())
    end
  end
  d.monitor = values
  d
end


#################### Logical Methods ####################

@promote_scalarvariate ScalarLogical


#################### Logical Constructors ####################

function Logical(expr::Expr, monitor::Union{Bool,Vector{Int}}=true)
  value = Float64(NaN)
  l = ScalarLogical(value, :nothing, 0, Int[], depfx(expr), depsrc(expr),
                    Symbol[])
  setmonitor!(l, monitor)
end

function Logical(d::Integer, expr::Expr, monitor::Union{Bool,Vector{Int}}=true)
  value = Array(Float64, fill(0, d)...)
  l = ArrayLogical(value, :nothing, 0, Int[], depfx(expr), depsrc(expr),
                   Symbol[])
  setmonitor!(l, monitor)
end


#################### Logical Updating ####################

function setinits!(l::AbstractLogical, m::Model, ::Any=nothing)
  l.value = l.eval(m)
  setmonitor!(l, l.monitor)
end

function update!(l::AbstractLogical, m::Model)
  l[:] = l.eval(m)
  l
end


#################### Stochastic Methods ####################

@promote_scalarvariate ScalarStochastic

function Base.showall(io::IO, s::AbstractStochastic)
  show(io, s)
  print(io, "\n\nDistribution:\n")
  show(io, s.distr)
  print(io, "\nFunction:\n")
  show(io, s.eval.code)
  print(io, "\n\nSource Nodes:\n")
  show(io, s.sources)
  print(io, "\n\nTarget Nodes:\n")
  show(io, s.targets)
end


#################### Stochastic Constructors ####################

function Stochastic(expr::Expr, monitor::Union{Bool,Vector{Int}}=true)
  value = Float64(NaN)
  s = ScalarStochastic(value, :nothing, 0, Int[], depfx(expr), depsrc(expr),
                       Symbol[], NullUnivariateDistribution())
  setmonitor!(s, monitor)
end

function Stochastic(d::Integer, expr::Expr,
                    monitor::Union{Bool,Vector{Int}}=true)
  value = Array(Float64, fill(0, d)...)
  s = ArrayStochastic(value, :nothing, 0, Int[], depfx(expr), depsrc(expr),
                      Symbol[], NullUnivariateDistribution())
  setmonitor!(s, monitor)
end


#################### Stochastic Updating ####################

function setinits!(s::ScalarStochastic, m::Model, x)
  s.value = x
  s.distr = s.eval(m)
  setmonitor!(s, s.monitor)
end

function setinits!(s::ArrayStochastic, m::Model, x)
  s.value = oftype(s.value, copy(x))
  s.distr = s.eval(m)
  if !isa(s.distr, UnivariateDistribution) && dims(s) != dims(s.distr)
    error("incompatible stochastic node and distribution structure dimensions")
  end
  setmonitor!(s, s.monitor)
end

function update!(s::AbstractStochastic, m::Model)
  s.distr = s.eval(m)
  s
end


#################### Stochastic Distribution Methods ####################

function invlink(s::AbstractStochastic, x, transform::Bool=true)
  invlink(s.distr, x, transform)
end

function link(s::AbstractStochastic, x, transform::Bool=true)
  link(s.distr, x, transform)
end

function logpdf(s::AbstractStochastic, transform::Bool=false)
  sum(logpdf(s.distr, s.value, transform))
end


#################### Auxiliary Functions ####################

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

function depfx(expr::Expr)
  eval(Expr(:function, :(model::Mamba.Model,), expr))
end
