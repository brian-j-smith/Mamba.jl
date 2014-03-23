#################### MCMCParam Methods ####################

function Base.show(io::IO, p::MCMCParam)
  msg = string(ifelse(p.monitor, "A ", "An un"),
               "monitored node of type \"", summary(p), "\"\n")
  print(io, msg)
  show(io, p.data)
  print(io, "\n")
end

function Base.showall(io::IO, p::MCMCParam)
  show(io, p)
  print(io, "\nFunction:\n")
  show(io, p.eval.code)
  print(io, "\n\nNode Dependencies:\n")
  show(io, p.deps)
  print(io, "\n")
end


#################### MCMCLogical Constructors ####################

function MCMCLogical(data, expr::Expr, monitor::Bool)
  MCMCLogical(data, monitor, paramfx(expr), paramdeps(expr))
end

function MCMCLogical(expr::Expr, monitor::Bool=true)
  data = convert(VariateType, NaN)
  MCMCLogical(data, expr, monitor)
end

function MCMCLogical(length::Integer, expr::Expr, monitor::Bool=true)
  data = Array(VariateType, length)
  fill!(data, NaN)
  MCMCLogical(data, expr, monitor)
end

function MCMCLogical(m::Integer, n::Integer, expr::Expr, monitor::Bool=true)
  data = Array(VariateType, m, n)
  fill!(data, NaN)
  MCMCLogical(data, expr, monitor)
end


#################### MCMCLogical Methods ####################

function initchain!(l::MCMCLogical, m::MCMCModel, chain::Integer)
  update!(l, m)
end

function logpdf(l::MCMCLogical)
  0
end

function update!(l::MCMCLogical, m::MCMCModel)
  l[:] = l.eval(m)
  l
end


#################### MCMCStochastic Constructors ####################

function MCMCStochastic{T}(data::T, expr::Expr, monitor::Bool)
  MCMCStochastic(data, monitor, paramfx(expr), paramdeps(expr), NullDistribution(), Array(T,0))
end

function MCMCStochastic(expr::Expr, monitor::Bool=true)
  data = convert(VariateType, NaN)
  MCMCStochastic(data, expr, monitor)
end

function MCMCStochastic(length::Integer, expr::Expr, monitor::Bool=true)
  data = Array(VariateType, length)
  fill!(data, NaN)
  MCMCStochastic(data, expr, monitor)
end

function MCMCStochastic(m::Integer, n::Integer, expr::Expr, monitor::Bool=true)
  data = Array(VariateType, m, n)
  fill!(data, NaN)
  MCMCStochastic(data, expr, monitor)
end


#################### MCMCStochastic Methods ####################

function Base.showall(io::IO, s::MCMCStochastic)
  show(io, s)
  print(io, "\n\nDistribution:\n")
  show(io, s.distr)
  print(io, "\nFunction:\n")
  show(io, s.eval.code)
  print(io, "\n\nNode Dependencies:\n")
  show(io, s.deps)
  print(io, "\n")
end

function initchain!(s::MCMCStochastic, m::MCMCModel, chain::Integer)
  length(s.inits) > 0 || error("missing initial values for stochastic node")
  i = (chain - 1) % length(s.inits) + 1
  s[:] = s.inits[i]
  update!(s, m)
  if isa(s.distr, Array) && size(s.data) != size(s.distr)
    error("stochastic parameter and distribution dimensions must match")
  end
  s
end

function insupport(s::MCMCStochastic)
  if isa(s.distr, Array)
    all(map(i -> insupport(s.distr[i], s.data[i]), 1:length(s)))
  else
    insupport(s.distr, s.data)
  end
end

function logpdf(s::MCMCStochastic)
  if isa(s.distr, Array)
    mapreduce(i -> logpdf(s.distr[i], s.data[i]), +, 1:length(s))
  else
    logpdf(s.distr, s.data)
  end
end

function update!(s::MCMCStochastic, m::MCMCModel)
  s.distr = s.eval(m)
  s
end


#################### Utility Functions ####################

function paramdeps(expr::Expr)
  if expr.head == :ref && expr.args[1] == :model && isa(expr.args[2], String)
    String[expr.args[2]]
  else
    mapreduce(paramdeps, union, expr.args)
  end
end

function paramdeps(expr)
  String[]
end

function paramfx(expr::Expr)
  eval(Main, Expr(:function, :(model::MCMCModel,), expr))
end
