#################### Sampler ####################

const samplerfxargs = [(:model, Mamba.Model), (:block, Integer)]


#################### Types and Constructors ####################

type NullFunction end

type SamplingBlock
  model::Model
  index::Int
  transform::Bool

  SamplingBlock(model::Model, index::Integer=0, transform::Bool=false) =
    new(model, index, transform)
end


Sampler(param::Symbol, args...) = Sampler([param], args...)

function Sampler(params::Vector{Symbol}, f::Function, tune::Any=Dict())
  Sampler(params, modelfx(samplerfxargs, f), tune, Symbol[])
end


function SamplerVariate{T<:SamplerTune, U<:Real}(x::AbstractVector{U}, tune::T)
  SamplerVariate{T}(x, tune)
end

function SamplerVariate(block::SamplingBlock, pargs...; kargs...)
  m = block.model
  SamplerVariate(unlist(block), m.samplers[block.index], m.iter, pargs...;
                 kargs...)
end

function SamplerVariate{T<:SamplerTune, U<:Real}(x::AbstractVector{U},
                                                 s::Sampler{T}, iter::Integer,
                                                 pargs...; kargs...)
  if iter == 1
    v = SamplerVariate{T}(x, pargs...; kargs...)
    s.tune = v.tune
  else
    v = SamplerVariate(x, s.tune)
  end
  v
end


#################### Base Methods ####################

function Base.show(io::IO, s::Sampler)
  print(io, "An object of type \"$(summary(s))\"\n")
  print(io, "Sampling Block Nodes:\n")
  show(io, s.params)
  print(io, "\n\n")
  show(io, "text/plain", first(code_typed(s.eval)))
  println(io)
end

function Base.showall(io::IO, s::Sampler)
  show(io, s)
  print(io, "\nTuning Parameters:\n")
  show(io, s.tune)
  print(io, "\n\nTarget Nodes:\n")
  show(io, s.targets)
end


#################### Variate Validators ####################

validate(v::SamplerVariate) = v

function validatebinary(v::SamplerVariate)
  all(insupport(Bernoulli, v)) ||
    throw(ArgumentError("variate is not a binary vector"))
  v
end

function validatesimplex(v::SamplerVariate)
  isprobvec(v) || throw(ArgumentError("variate is not a probability vector"))
  v
end


#################### sample! Generics ####################

function sample!(v::SamplerVariate, density::Nullable; args...)
  isnull(density) && error("must specify a target density in $(typeof(v))",
                           " constructor or sample! method")
  sample!(v, get(density); args...)
end


#################### Simulation Methods ####################

function gradlogpdf!{T<:Real}(block::SamplingBlock, x::AbstractArray{T},
                              dtype::Symbol=:forward)
  gradlogpdf!(block.model, x, block.index, block.transform, dtype=dtype)
end

function logpdf!{T<:Real}(block::SamplingBlock, x::AbstractArray{T})
  logpdf!(block.model, x, block.index, block.transform)
end

function logpdfgrad!{T<:Real}(block::SamplingBlock, x::AbstractVector{T},
                              dtype::Symbol)
  grad = gradlogpdf!(block, x, dtype)
  logf = logpdf!(block, x)
  (logf, ifelse(isfinite(grad), grad, 0.0))
end

function unlist(block::SamplingBlock)
  unlist(block.model, block.index, block.transform)
end

function relist{T<:Real}(block::SamplingBlock, x::AbstractArray{T})
  relist(block.model, x, block.index, block.transform)
end


#################### Auxiliary Functions ####################

asvec(x::Union{Number, Symbol}) = [x]
asvec(x::Vector) = x
