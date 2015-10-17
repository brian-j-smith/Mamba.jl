#################### Dependent ####################

#################### Base Methods ####################

function Base.show(io::IO, d::AbstractDependent)
  msg = string(ifelse(isempty(d.monitor), "An un", "A "),
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

function names(d::AbstractDependent)
  names(d, d.symbol)
end

function setmonitor!(d::AbstractDependent, monitor::Bool)
  value = monitor ? Int[0] : Int[]
  setmonitor!(d, value)
end

function setmonitor!(d::AbstractDependent, monitor::Vector{Int})
  values = monitor
  n = length(unlist(d))
  if n > 0 && !isempty(monitor)
    if monitor[1] == 0
      values = collect(1:n)
    elseif minimum(monitor) < 1 || maximum(monitor) > n
      throw(BoundsError())
    end
  end
  d.monitor = values
  d
end


#################### Distribution Fallbacks ####################

unlist(d::AbstractDependent, transform::Bool=false) =
  unlist(d, d.value, transform)

unlist(d::AbstractDependent, x::Real, transform::Bool=false) = [x]

unlist(d::AbstractDependent, x::AbstractArray, transform::Bool=false) = vec(x)

relist(d::AbstractDependent, x::AbstractArray, transform::Bool=false) =
  relistlength(d, x, transform)[1]

logpdf(d::AbstractDependent, transform::Bool=false) = 0.0

logpdf(d::AbstractDependent, x, transform::Bool=false) = 0.0


#################### Logical ####################

@promote_scalarvariate ScalarLogical


#################### Constructors ####################

function Logical(expr::Expr, monitor::Union{Bool,Vector{Int}}=true)
  value = Float64(NaN)
  l = ScalarLogical(value, :nothing, Int[], depfx(expr), depsrc(expr), Symbol[])
  setmonitor!(l, monitor)
end

function Logical(d::Integer, expr::Expr, monitor::Union{Bool,Vector{Int}}=true)
  value = Array(Float64, fill(0, d)...)
  l = ArrayLogical(value, :nothing, Int[], depfx(expr), depsrc(expr), Symbol[])
  setmonitor!(l, monitor)
end

function Logical(f::Function, monitor::Union{Bool,Vector{Int}}=true)
  Logical(modelexpr(f), monitor)
end

function Logical(d::Integer, f::Function, monitor::Union{Bool,Vector{Int}}=true)
  Logical(d, modelexpr(f), monitor)
end


#################### Updating ####################

function setinits!(l::AbstractLogical, m::Model, ::Any=nothing)
  l.value = l.eval(m)
  setmonitor!(l, l.monitor)
end

function update!(l::AbstractLogical, m::Model)
  l[:] = l.eval(m)
  l
end


#################### Distribution Methods ####################

relistlength(d::ScalarLogical, x::AbstractArray, transform::Bool=false) =
  (x[1], 1)

function relistlength(d::ArrayLogical, x::AbstractArray, transform::Bool=false)
  n = length(d)
  value = reshape(x[1:n], size(d))
  (value, n)
end


#################### Stochastic ####################

#################### Base Methods ####################

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


#################### Constructors ####################

function Stochastic(expr::Expr, monitor::Union{Bool,Vector{Int}}=true)
  value = Float64(NaN)
  s = ScalarStochastic(value, :nothing, Int[], depfx(expr), depsrc(expr),
                       Symbol[], NullUnivariateDistribution())
  setmonitor!(s, monitor)
end

function Stochastic(d::Integer, expr::Expr,
                    monitor::Union{Bool,Vector{Int}}=true)
  value = Array(Float64, fill(0, d)...)
  s = ArrayStochastic(value, :nothing, Int[], depfx(expr), depsrc(expr),
                      Symbol[], NullUnivariateDistribution())
  setmonitor!(s, monitor)
end

function Stochastic(f::Function, monitor::Union{Bool,Vector{Int}}=true)
  Stochastic(modelexpr(f), monitor)
end

function Stochastic(d::Integer, f::Function,
                    monitor::Union{Bool,Vector{Int}}=true)
  Stochastic(d, modelexpr(f), monitor)
end


#################### Updating ####################

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


#################### Distribution Methods ####################

function unlist(s::AbstractStochastic, transform::Bool=false)
  unlist(s, s.value, transform)
end

function unlist(s::AbstractStochastic, x::Real, transform::Bool=false)
  unlist(s, [x], transform)
end

function unlist(s::AbstractStochastic, x::AbstractArray, transform::Bool=false)
  transform ? unlist_sub(s.distr, link_sub(s.distr, x)) :
              unlist_sub(s.distr, x)
end

function relist(s::AbstractStochastic, x::AbstractArray, transform::Bool=false)
  relistlength(s, x, transform)[1]
end

function relistlength(s::AbstractStochastic, x::AbstractArray,
                      transform::Bool=false)
  value, n = relistlength_sub(s.distr, s, x)
  (transform ? invlink_sub(s.distr, value) : value, n)
end

function logpdf(s::AbstractStochastic, transform::Bool=false)
  logpdf(s, s.value, transform)
end

function logpdf(s::AbstractStochastic, x, transform::Bool=false)
  logpdf_sub(s.distr, x, transform)
end

rand(s::AbstractStochastic) = rand_sub(s.distr, s.value)


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
