#################### Binary MCMC Model Composition ####################

#################### Types and Constructors ####################

const BMC3Form = Union{Int, Vector{Vector{Int}}}

mutable struct BMC3Tune{F<:BMC3Form} <: SamplerTune
  logf::Union{Function, Missing}
  k::F

  BMC3Tune{F}() where {F<:BMC3Form} = new{F}()

  BMC3Tune{F}(x::Vector, k::F) where {F<:BMC3Form} = new{F}(missing, k)

  BMC3Tune{F}(x::Vector, k::F, logf::Function) where {F<:BMC3Form} =
    new{F}(logf, k)
end


const BMC3IntVariate = SamplerVariate{BMC3Tune{Int}}
const BMC3VecVariate = SamplerVariate{BMC3Tune{Vector{Vector{Int}}}}

BMC3Variate(x::Vector, logf::Function; k::F=1) where {F<:BMC3Form} =
  SamplerVariate{BMC3Tune{F}}(x, k, logf)


function validate(v::BMC3IntVariate)
  n = length(v)
  v.tune.k <= n || throw(ArgumentError("k exceeds variate length $n"))
  validatebinary(v)
end

function validate(v::BMC3VecVariate)
  n = length(v)
  mapreduce(maximum, max, v.tune.k) <= n ||
    throw(ArgumentError("indices in k exceed variate length $n"))
  validatebinary(v)
end


#################### Sampler Constructor ####################

function BMC3(params::ElementOrVector{Symbol}; k::F=1) where {F<:BMC3Form}
  samplerfx = function(model::Model, block::Integer)
    block = SamplingBlock(model, block)
    v = SamplerVariate(block, k)
    sample!(v, x -> logpdf!(block, x))
    relist(block, v)
  end
  Sampler(params, samplerfx, BMC3Tune{F}())
end


#################### Sampling Functions ####################

sample!(v::SamplerVariate{BMC3Tune{F}}) where {F<:BMC3Form} = sample!(v, v.tune.logf)

function sample!(v::SamplerVariate{BMC3Tune{F}}, logf::Function) where {F<:BMC3Form}
  x = v[:]
  idx = randind(v)
  x[idx] = 1.0 - v[idx]
  if rand() < exp(logf(x) - logf(v.value))
    v[:] = x
  end
  v
end

randind(v::BMC3IntVariate) = sample(1:length(v), v.tune.k, replace=false)
randind(v::BMC3VecVariate) = sample(v.tune.k)
