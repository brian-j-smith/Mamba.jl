#################### Binary MCMC Model Composition ####################

#################### Types and Constructors ####################

type BMC3Tune <: SamplerTune
  logf::Nullable{Function}
  k::Int
  indexset::Vector{Vector{Int}}

  BMC3Tune() = new()

  BMC3Tune(x::Vector, logf::Nullable{Function}; k::Integer=1) =
    new(logf, k, Vector{Vector{Int}}())
end

BMC3Tune(x::Vector; args...) =
   BMC3Tune(x, Nullable{Function}(); args...)

BMC3Tune(x::Vector, logf::Function; args...) =
   BMC3Tune(x, Nullable{Function}(logf); args...)


typealias BMC3Variate SamplerVariate{BMC3Tune}

function validate(v::BMC3Variate)
  n = length(v)
  v.tune.k <= n || throw(ArgumentError("k exceeds variate length $n"))
  validatebinary(v)
end


#################### Sampler Constructor ####################

function BMC3(params::ElementOrVector{Symbol}; args...)
  samplerfx = function(model::Model, block::Integer)
    block = SamplingBlock(model, block)
    v = SamplerVariate(block; args...)
    sample!(v, x -> logpdf!(block, x))
    relist(block, v)
  end
  Sampler(params, samplerfx, BMC3Tune())
end


#################### Sampling Functions ####################

sample!(v::BMC3Variate) = sample!(v, v.tune.logf)

function sample!(v::BMC3Variate, logf::Function)
  x = v[:]
  idx = randperm(length(v))[1:v.tune.k]
  x[idx] = 1.0 - v[idx]
  if rand() < exp(logf(x) - logf(v.value))
    v[:] = x
  end
  v
end


#################### Legacy Sampler Code ####################

function BMC3(params::ElementOrVector{Symbol}, indexset::Vector{Vector{Int}})
  samplerfx = function(model::Model, block::Integer)
    block = SamplingBlock(model, block)
    v = SamplerVariate(block, 0)
    bmc3!(v, indexset, x -> logpdf!(block, x))
    relist(block, v)
  end
  Sampler(params, samplerfx, BMC3Tune())
end


function bmc3!(v::BMC3Variate, logf::Function; k::Integer=1)
  v.tune.k = k
  sample!(v, logf)
end

function bmc3!(v::BMC3Variate, indexset::Vector{Vector{Int}}, logf::Function)
  v.tune.indexset = indexset
  x = v[:]
  idx = rand(v.tune.indexset)
  x[idx] = 1.0 - v[idx]
  if rand() < exp(logf(x) - logf(v.value))
    v[:] = x
  end
  v
end
