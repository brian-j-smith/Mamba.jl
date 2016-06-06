############ Binary Metropolised Gibbs Sampler ##############

#################### Types ####################

type BMGTune <: SamplerTune
  logf::Nullable{Function}
  k::Int

  BMGTune() = new()

  BMGTune(x::Vector, logf::Nullable{Function}; k::Integer=1) = new(logf, k)
end

BMGTune(x::Vector; args...) =
  BMGTune(x, Nullable{Function}(); args...)

BMGTune(x::Vector, logf::Function; args...) =
  BMGTune(x, Nullable{Function}(logf); args...)


typealias BMGVariate SamplerVariate{BMGTune}

function validate(v::BMGVariate)
  n = length(v)
  v.tune.k <= n || throw(ArgumentError("k exceeds variate length $n"))
  validatebinary(v)
end


#################### Sampler Constructor ####################

function BMG(params::ElementOrVector{Symbol}; args...)
  samplerfx = function(model::Model, block::Integer)
    block = SamplingBlock(model, block)
    v = SamplerVariate(block; args...)
    sample!(v, x -> logpdf!(block, x))
    relist(block, v)
  end
  Sampler(params, samplerfx, BMGTune())
end


#################### Sampling Functions ####################

sample!(v::BMGVariate) = sample!(v, v.tune.logf)

function sample!(v::BMGVariate, logf::Function)
  n = length(v)
  probs = Vector{Float64}(n)
  idx = sample(1:n, v.tune.k, replace=false)

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
