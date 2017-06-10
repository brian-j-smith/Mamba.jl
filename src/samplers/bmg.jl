#################### Binary Metropolised Gibbs Sampler ####################

#################### Types and Constructors ####################

const BMGForm = Union{Int, Vector{Vector{Int}}}

type BMGTune{F<:BMGForm} <: SamplerTune
  logf::Nullable{Function}
  k::F

  BMGTune() = new()

  BMGTune(x::Vector, k::F) = new(Nullable{Function}(), k)

  BMGTune(x::Vector, k::F, logf::Function) = new(Nullable{Function}(logf), k)
end


const BMGIntVariate = SamplerVariate{BMGTune{Int}}
const BMGVecVariate = SamplerVariate{BMGTune{Vector{Vector{Int}}}}

BMGVariate{F<:BMGForm}(x::Vector, logf::Function; k::F=1) =
  SamplerVariate{BMGTune{F}}(x, k, logf)


function validate(v::BMGIntVariate)
  n = length(v)
  v.tune.k <= n || throw(ArgumentError("k exceeds variate length $n"))
  validatebinary(v)
end

function validate(v::BMGVecVariate)
  n = length(v)
  mapreduce(maximum, max, v.tune.k) <= n ||
    throw(ArgumentError("indices in k exceed variate length $n"))
  validatebinary(v)
end


#################### Sampler Constructor ####################

function BMG{F<:BMGForm}(params::ElementOrVector{Symbol}; k::F=1)
  samplerfx = function(model::Model, block::Integer)
    block = SamplingBlock(model, block)
    v = SamplerVariate(block, k)
    sample!(v, x -> logpdf!(block, x))
    relist(block, v)
  end
  Sampler(params, samplerfx, BMGTune{F}())
end


#################### Sampling Functions ####################

sample!{F<:BMGForm}(v::SamplerVariate{BMGTune{F}}) = sample!(v, v.tune.logf)

function sample!{F<:BMGForm}(v::SamplerVariate{BMGTune{F}}, logf::Function)
  n = length(v)
  probs = Vector{Float64}(n)
  idx = randind(v)

  pbernoulli! = function(x, probs)
    for i in idx
      x_i = x[i]

      x[i] = 0.0
      logf0 = logf(x)
      x[i] = 1.0
      logf1 = logf(x)

      p = invlogit(logf1 - logf0)
      probs[i] = 0.0 < p < 1.0 ? p : 0.5

      x[i] = x_i
    end
    probs
  end

  x = v[:]
  pbernoulli!(x, probs)
  theta = map(p -> rand() < p, probs[idx])

  if n == 1
    v[idx] = theta
  else
    y = v[:]
    y[idx] = theta
    qy = mapreduce(i -> y[i] == 1.0 ? log(probs[i]) : log(1.0 - probs[i]), +,
                   idx)

    pbernoulli!(y, probs)
    qx = mapreduce(i -> x[i] == 1.0 ? log(probs[i]) : log(1.0 - probs[i]), +,
                   idx)

    if rand() < exp((logf(y) - qy) - (logf(x) - qx))
      v[idx] = theta
    end
  end

  v
end

randind(v::BMGIntVariate) = sample(1:length(v), v.tune.k, replace=false)
randind(v::BMGVecVariate) = sample(v.tune.k)
