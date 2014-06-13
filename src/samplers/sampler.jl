#################### MCMCSampler Constructor ####################

function MCMCSampler(params::Vector{Symbol}, expr::Expr, tune::Dict=Dict())
  MCMCSampler(params, samplerfx(expr), tune, Symbol[])
end


#################### MCMCSampler Methods ####################

function Base.show(io::IO, s::MCMCSampler)
  print(io, "An object of type \"$(summary(s))\"\n")
  print(io, "Sampling Block Nodes:\n")
  show(io, s.params)
  print(io, "\n\n")
  show(io, s.eval.code)
  println(io)
end

function Base.showall(io::IO, s::MCMCSampler)
  show(io, s)
  print(io, "\nTuning Parameters:\n")
  show(io, s.tune)
  print(io, "\n\nTarget Nodes:\n")
  show(io, s.targets)
end


#################### Utility Functions ####################

function samplerfx(expr::Expr)
  eval(Main, Expr(:function, :(model::MCMCModel, block::Integer), expr))
end
