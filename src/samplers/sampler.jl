#################### Sampler ####################

const samplerfxargs = [(:model, :Model), (:block, :Integer)]


#################### Constructors ####################

function Sampler(params::Vector{Symbol}, f::Function, tune::Dict=Dict())
  Sampler(params, modelfx(samplerfxargs, f), tune, Symbol[])
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
