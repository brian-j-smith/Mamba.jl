############ Binary Metropolised Gibbs Sampler ##############

#################### Types ####################

type BMGTune <: SamplerTune
  BMGTune(d::Integer=0) = new()
end

type BMGVariate <: SamplerVariate
  value::Vector{Float64}
  tune::BMGTune

  function BMGVariate{T<:Real}(x::AbstractVector{T}, tune::BMGTune)
    all(insupport(Bernoulli, x)) ||
      throw(ArgumentError("x is not a binary vector"))
    new(x, tune)
  end
end

function BMGVariate{T<:Real}(x::AbstractVector{T}, tune=nothing)
  BMGVariate(x, BMGTune(length(x)))
end


#################### Sampler Constructor ####################

function BMG(params::Vector{Symbol}; k::Integer=1)
  samplerfx = function(model::Model, block::Integer)
    v = variate!(BMGVariate, unlist(model, block),
                 model.samplers[block], model.iter)
    f = x -> logpdf!(model, x, block)
    bmg!(v, f, k=k)
    relist(model, v, block)
  end
  Sampler(params, samplerfx, BMGTune())
end


#################### Sampling Functions ####################

function bmg!(v::BMGVariate, logf::Function; k::Integer=1)
  d = length(v)
  k <= d || throw(ArgumentError("k is greater than length(v)"))

  probs = Array(Float64, d)
  idx = randperm(d)[1:k]

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

  if d == 1
    v[idx] = theta
  else
    y = v[:]
    y[idx] = theta
    qy = mapreduce(i -> y[i] == 1.0 ? log(probs[i]) : log(1.0 - probs[i]), +, idx)

    pbernoulli!(y, probs)
    qx = mapreduce(i -> x[i] == 1.0 ? log(probs[i]) : log(1.0 - probs[i]), +, idx)

    if rand() < exp((logf(y) - qy) - (logf(x) - qx))
      v[idx] = theta
    end
  end

  v
end
