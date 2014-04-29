#################### MCMCDependent Methods ####################

function Base.show(io::IO, d::MCMCDependent)
  msg = string(ifelse(any(d.monitor), "A ", "An un"),
               "monitored node of type \"", summary(d), "\"\n")
  print(io, msg)
  show(io, d.value)
  print(io, "\n")
end

function Base.showall(io::IO, d::MCMCDependent)
  show(io, d)
  print(io, "\nFunction:\n")
  show(io, d.eval.code)
  print(io, "\n\nNode Dependencies:\n")
  show(io, d.sources)
  print(io, "\n")
end

identity(d::MCMCDependent, x) = x

invlink(d::MCMCDependent, x) = x

link(d::MCMCDependent, x) = x

logpdf(d::MCMCDependent, transform::Bool=false) = 0.0

function setmonitor!(d::MCMCDependent, monitor::Union(Bool,Vector{Bool}))
  if isa(monitor, Bool)
    values = fill(monitor, length(d))
  elseif length(monitor) == length(d)
    values = deepcopy(monitor)
  else
    error("node and monitor dimensions must match")
  end
  d.monitor = values
  d
end


#################### MCMCLogical Constructors ####################

function MCMCLogical(value, expr::Expr, monitor::Union(Bool,Vector{Bool}))
  l = MCMCLogical(value, String[], Bool[], paramfx(expr), paramsrc(expr))
  setmonitor!(l, monitor)
end

function MCMCLogical(expr::Expr, monitor::Union(Bool,Vector{Bool})=true)
  value = convert(VariateType, NaN)
  MCMCLogical(value, expr, monitor)
end

function MCMCLogical(length::Integer, expr::Expr,
           monitor::Union(Bool,Vector{Bool})=true)
  value = Array(VariateType, length)
  fill!(value, NaN)
  MCMCLogical(value, expr, monitor)
end

function MCMCLogical(m::Integer, n::Integer, expr::Expr,
           monitor::Union(Bool,Vector{Bool})=true)
  value = Array(VariateType, m, n)
  fill!(value, NaN)
  MCMCLogical(value, expr, monitor)
end


#################### MCMCLogical Methods ####################

setinits!(l::MCMCLogical, m::MCMCModel, x=nothing) = update!(l, m)

function update!(l::MCMCLogical, m::MCMCModel)
  l[:] = l.eval(m)
  l
end


#################### MCMCStochastic Constructors ####################

function MCMCStochastic{T}(value::T, expr::Expr,
           monitor::Union(Bool,Vector{Bool}))
  s = MCMCStochastic(value, String[], Bool[], paramfx(expr), paramsrc(expr),
                     NullDistribution())
  setmonitor!(s, monitor)
end

function MCMCStochastic(expr::Expr, monitor::Union(Bool,Vector{Bool})=true)
  value = convert(VariateType, NaN)
  MCMCStochastic(value, expr, monitor)
end

function MCMCStochastic(length::Integer, expr::Expr,
           monitor::Union(Bool,Vector{Bool})=true)
  value = Array(VariateType, length)
  fill!(value, NaN)
  MCMCStochastic(value, expr, monitor)
end

function MCMCStochastic(m::Integer, n::Integer, expr::Expr,
           monitor::Union(Bool,Vector{Bool})=true)
  value = Array(VariateType, m, n)
  fill!(value, NaN)
  MCMCStochastic(value, expr, monitor)
end


#################### MCMCStochastic Methods ####################

function Base.showall(io::IO, s::MCMCStochastic)
  show(io, s)
  print(io, "\n\nDistribution:\n")
  show(io, s.distr)
  print(io, "\nFunction:\n")
  show(io, s.eval.code)
  print(io, "\n\nNode Dependencies:\n")
  show(io, s.sources)
  print(io, "\n")
end

function setinits!(s::MCMCStochastic, m::MCMCModel, x)
  s[:] = convert(typeof(s.value), x)
  update!(s, m)
  if isa(s.distr, Array) && size(s.value) != size(s.distr)
    error("stochastic parameter and distribution dimensions must match")
  end
  s
end

insupport(s::MCMCStochastic) = all(mapdistr(insupport, s, s.value))

invlink(s::MCMCStochastic, x) = mapdistr(invlink, s, x)

link(s::MCMCStochastic, x) =  mapdistr(link, s, x)

function logpdf(s::MCMCStochastic, transform::Bool=false)
  sum(mapdistr(logpdf, s, s.value, transform))
end

function mapdistr(f::Function, s::MCMCStochastic, x, args...)
  if isa(s.distr, Array)
    n = length(s.distr)
    length(x) == n ||
      error("dimension of stochastic distribution and argument must match")
    map(i -> f(s.distr[i], x[i], args...), 1:n)
  else
    f(s.distr, x, args...)
  end
end

function update!(s::MCMCStochastic, m::MCMCModel)
  s.distr = s.eval(m)
  s
end


#################### Utility Functions ####################

function paramsrc(expr::Expr)
  if expr.head == :ref && expr.args[1] == :model && isa(expr.args[2], String)
    String[expr.args[2]]
  else
    mapreduce(paramsrc, union, expr.args)
  end
end

function paramsrc(expr)
  String[]
end

function paramfx(expr::Expr)
  eval(Main, Expr(:function, :(model::MCMCModel,), expr))
end
