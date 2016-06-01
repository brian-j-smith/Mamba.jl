#################### Binary MCMC Model Composition ####################

#################### Types and Constructors ####################

typealias BMC3Form Union{Int, Vector{Vector{Int}}}

type BMC3Tune{F<:BMC3Form} <: SamplerTune
  logf::Nullable{Function}
  indices::F
  n::Int
  m::Vector{Float64}
  v::Vector{Float64}
  epsilon::Float64

  BMC3Tune() = new()

  BMC3Tune(x::Vector, indices::F; args...) =
    BMC3Tune{F}(x, indices, Nullable{Function}(); args...)

  BMC3Tune(x::Vector, indices::F, logf::Function; args...) =
    BMC3Tune{F}(x, indices, Nullable{Function}(logf); args...)

  BMC3Tune(x::Vector, indices::F, logf::Nullable{Function};
           epsilon::Real=1 / length(x)) =
    new(logf, indices, 0, zeros(x), zeros(x), epsilon)
end


BMC3Variate(x::Vector; k::Integer=1, args...) =
   SamplerVariate{BMC3Tune{Int}}(x, k; args...)

BMC3Variate(x::Vector, logf::Function; k::Integer=1, args...) =
   SamplerVariate{BMC3Tune{Int}}(x, k, logf; args...)

BMC3Variate(x::Vector, indexset::Vector{Vector{Int}}; args...) =
  SamplerVariate{BMC3Tune{Vector{Vector{Int}}}}(x, indexset; args...)

BMC3Variate(x::Vector, indexset::Vector{Vector{Int}}, logf::Function; args...) =
  SamplerVariate{BMC3Tune{Vector{Vector{Int}}}}(x, indexset, logf; args...)


function validate(v::SamplerVariate{BMC3Tune{Int}})
  n = length(v)
  v.tune.indices <= n ||
    throw(ArgumentError("indices exceed variate length $n"))
  validatebinary(v)
end

function validate(v::SamplerVariate{BMC3Tune{Vector{Vector{Int}}}})
  n = length(v)
  mapreduce(maximum, max, v.tune.indices) <= n ||
    throw(ArgumentError("indices exceed variate length $n"))
  validatebinary(v)
end


#################### Sampler Constructor ####################

function BMC3(params::ElementOrVector{Symbol}; k::Integer=1, args...)
  BMC3Sampler(params, k; args...)
end

function BMC3(params::ElementOrVector{Symbol}, indexset::Vector{Vector{Int}};
              args...)
  BMC3Sampler(params, indexset; args...)
end

function BMC3Sampler{F<:BMC3Form}(params::ElementOrVector{Symbol}, indices::F;
                                  adapt::Symbol=:none, atype::Symbol=:var,
                                  args...)
  adapt in [:all, :burnin, :none] ||
    throw(ArgumentError("adapt must be one of :all, :burnin, or :none"))
  atype in [:var, :freq] ||
    throw(ArgumentError("atype must be one of :var or :freq"))

  samplerfx = function(model::Model, block::Integer)
    block = SamplingBlock(model, block)
    v = SamplerVariate(block, indices; args...)
    isadapt = adapt == :burnin ? model.iter <= model.burnin :
              adapt == :all ? true : false
    sample!(v, x -> logpdf!(block, x), adapt=isadapt, atype=atype)
    relist(block, v)
  end
  Sampler(params, samplerfx, BMC3Tune{F}())
end


#################### Sampling Functions ####################

sample!{F<:BMC3Form}(v::SamplerVariate{BMC3Tune{F}}; args...) =
  sample!(v, v.tune.logf; args...)

function sample!{F<:BMC3Form}(v::SamplerVariate{BMC3Tune{F}}, logf::Function;
                              adapt::Bool=false, atype::Symbol=:var)
  x = v[:]
  setadapt!(v, adapt)
  idx = randind(v, atype)
  x[idx] = 1.0 - v[idx]
  if rand() < exp(logf(x) - logf(v.value))
    v[:] = x
  end
  v
end

function randind(v::SamplerVariate{BMC3Tune{Int}}, atype)
  wv =  atype == :var ?
        (1 - v.tune.epsilon) * (v.tune.v / (v.tune.n - 1)) + v.tune.epsilon :
        (1 - v.tune.epsilon) * v.tune.m + v.tune.epsilon
  sample(1:length(v), wv, v.tune.indices)
end

function randind(v::SamplerVariate{BMC3Tune{Vector{Vector{Int}}}})
  wv =  atype == :var ?
        (1 - v.tune.epsilon) * (v.tune.v / (v.tune.n - 1)) + v.tune.epsilon :
        (1 - v.tune.epsilon) * v.tune.m + v.tune.epsilon
  wvb = zeros(length(v.tune.indices))
  for i in 1:length(v.tune.indices)
    wvb[i] = mean(wv[v.tune.indices[i]])
  end
  sample(v.tune.indices, wvb)
end

function setadapt!(v::BMC3Variate, adapt::Bool)
  if adapt
    v.tune.n += 1
    d = v.value - v.tune.m
    v.tune.m += d / v.tune.n
    v.tune.v += d .* (x - v.tune.m)
  end
end

