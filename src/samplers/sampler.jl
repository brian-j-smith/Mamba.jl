#################### MCMCSampler Constructor ####################

function MCMCSampler{T<:String}(params::Vector{T}, expr::Expr,
           tune::Dict=Dict())
  MCMCSampler(String[params...], samplerfx(expr), tune, String[], String[])
end


#################### MCMCSampler Methods ####################

function Base.show(io::IO, s::MCMCSampler)
  print(io, "An object of type \"$(summary(s))\"\n")
  print(io, "Sampling Block Nodes:\n")
  show(io, s.params)
  print(io, "\n\n")
  show(io, s.eval.code)
  print(io, "\n")
end

function Base.showall(io::IO, s::MCMCSampler)
  show(io, s)
  print(io, "\nTuning Parameters:\n")
  show(io, s.tune)
  print(io, "\nSource Nodes:\n")
  show(io, s.sources)
  print(io, "\nTarget Nodes:\n")
  show(io, s.targets)
  print(io, "\n")
end


#################### Utility Functions ####################

function samplerfx(expr::Expr)
  eval(Main, Expr(:function, :(model::MCMCModel, block::Integer), expr))
end
