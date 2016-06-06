#################### Binary MCMC Model Composition ####################

#################### Types and Constructors ####################

typealias BMC3Form Union{Int, Vector{Vector{Int}}}

type BMC3Tune{F<:BMC3Form} <: SamplerTune
  logf::Nullable{Function}
  k::F

  BMC3Tune() = new()

  BMC3Tune(x::Vector, k::F) = new(Nullable{Function}(), k)

  BMC3Tune(x::Vector, k::F, logf::Function) = new(Nullable{Function}(logf), k)
end


typealias BMC3IntVariate SamplerVariate{BMC3Tune{Int}}
typealias BMC3VecVariate SamplerVariate{BMC3Tune{Vector{Vector{Int}}}}

BMC3Variate{F<:BMC3Form}(x::Vector, logf::Function; k::F=1) =
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

function BMC3{F<:BMC3Form}(params::ElementOrVector{Symbol}; k::F=1)
  samplerfx = function(model::Model, block::Integer)
    block = SamplingBlock(model, block)
    v = SamplerVariate(block, k)
    sample!(v, x -> logpdf!(block, x))
    relist(block, v)
  end
  Sampler(params, samplerfx, BMC3Tune{F}())
end


#################### Sampling Functions ####################

sample!{F<:BMC3Form}(v::SamplerVariate{BMC3Tune{F}}) = sample!(v, v.tune.logf)

function sample!{F<:BMC3Form}(v::SamplerVariate{BMC3Tune{F}}, logf::Function)
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
