#################### Sampler ####################

const samplerfxargs = [(:model, :Model), (:block, :Integer)]


#################### Constructors ####################

function SamplerVariate{T<:SamplerTune, U<:Real}(x::AbstractVector{U}, tune::T)
  SamplerVariate{T}(x, tune)
end

function SamplerVariate{T<:SamplerTune, U<:Real}(x::AbstractVector{U},
                                                 s::Sampler{T}, iter::Integer)
  if iter == 1
    v = SamplerVariate{T}(x)
    s.tune = v.tune
  else
    v = SamplerVariate{T}(x, s.tune)
  end
  v
end

function SamplerVariate(m::Model, block::Integer, transform::Bool=false)
  SamplerVariate(unlist(m, block, transform), m.samplers[block], m.iter)
end


function Sampler(params::Vector{Symbol}, f::Function, tune::Any=Dict())
  Sampler(params, modelfx(samplerfxargs, f), tune, Symbol[])
end


#################### Variate Validators ####################

validate(v::SamplerVariate) = v

function validate(v::Union{BHMCVariate, BMC3Variate, BMGVariate})
  all(insupport(Bernoulli, v)) ||
    throw(ArgumentError("variate is not a binary vector"))
  v
end

function validate(v::SliceSimplexVariate)
  isprobvec(v) || throw(ArgumentError("variate is not a probability vector"))
  v
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

function logpdfgrad{T<:Real}(m::Model, x::Vector{T}, block::Integer,
                             dtype::Symbol)
  logf = logpdf(m, x, block, true)
  grad = isfinite(logf) ?
           gradlogpdf(m, x, block, true, dtype=dtype) :
           zeros(x)
  logf, grad
end

function logpdfgrad!{T<:Real}(m::Model, x::Vector{T}, block::Integer,
                              dtype::Symbol)
  logf = logpdf!(m, x, block, true)
  grad = isfinite(logf) ?
           gradlogpdf!(m, x, block, true, dtype=dtype) :
           zeros(length(x))
  logf, grad
end
