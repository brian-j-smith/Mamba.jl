#################### Binary MCMC Model Composition ####################

#################### Types and Constructors ####################

type BMC3Tune <: SamplerTune
  indexset::Vector{Vector{Int}}

  function BMC3Tune(value::Vector{Float64}=Float64[])
    new(
      Vector{Vector{Int}}[]
    )
  end
end


typealias BMC3Variate SamplerVariate{BMC3Tune}


#################### Sampler Constructor ####################

function BMC3(params::ElementOrVector{Symbol}; k::Integer=1)
  samplerfx = function(model::Model, block::Integer)
    v = SamplerVariate(model, block)
    f = x -> logpdf!(model, x, block)
    bmc3!(v, f, k=k)
    relist(model, v, block)
  end
  Sampler(params, samplerfx, BMC3Tune())
end

function BMC3(params::ElementOrVector{Symbol}, indexset::Vector{Vector{Int}})
  samplerfx = function(model::Model, block::Integer)
    v = SamplerVariate(model, block)
    f = x -> logpdf!(model, x, block)
    bmc3!(v, indexset, f)
    relist(model, v, block)
  end
  Sampler(params, samplerfx, BMC3Tune())
end


#################### Sampling Functions ####################

function bmc3!(v::BMC3Variate, logf::Function; k::Integer=1)
  d = length(v)
  k <= d || throw(ArgumentError("k is greater than length(v)"))
  bmc3_sub!(v, randperm(d)[1:k], logf)
end

function bmc3!(v::BMC3Variate, indexset::Vector{Vector{Int}}, logf::Function)
  v.tune.indexset = indexset
  bmc3_sub!(v, indexset[rand(1:length(indexset))], logf)
end

function bmc3_sub!(v::BMC3Variate, idx::Vector{Int}, logf::Function)
  x = v[:]
  x[idx] = 1.0 - v[idx]
  if rand() < exp(logf(x) - logf(v.value))
    v[:] = x
  end
  v
end
