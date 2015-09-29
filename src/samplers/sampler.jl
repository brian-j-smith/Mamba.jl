#################### Sampler ####################

#################### Constructors ####################

function Sampler(params::Vector{Symbol}, expr::Expr, tune::Dict=Dict())
  Sampler(params, samplerfx(expr), tune, Symbol[])
end


#################### Base Methods ####################

function Base.show(io::IO, s::Sampler)
  print(io, "An object of type \"$(summary(s))\"\n")
  print(io, "Sampling Block Nodes:\n")
  show(io, s.params)
  print(io, "\n\n")
  show(io, s.eval.code)
  println(io)
end

function Base.showall(io::IO, s::Sampler)
  show(io, s)
  print(io, "\nTuning Parameters:\n")
  show(io, s.tune)
  print(io, "\n\nTarget Nodes:\n")
  show(io, s.targets)
end


#################### Auxiliary Functions ####################

function samplerfx(expr::Expr)
  eval(Expr(:function, :(model::Mamba.Model, block::Integer), expr))
end
