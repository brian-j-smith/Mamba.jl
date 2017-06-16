##################### Binary Individual Adaptation ####################

#################### Types and Constructors ####################

type BIATune <: SamplerTune
  logf::Nullable{Function}
  A::Vector{Float64}
  D::Vector{Float64}

  epsilon::Float64
  decay::Float64
  target::Float64

  iter::Int

  BIATune() = new()

  function BIATune(x::Vector, logf::Nullable{Function};
                   A::Vector{Float64} = ones(x) / length(x),
                   D::Vector{Float64} = ones(x) / length(x),
                   epsilon::Real = 0.01 / length(x), decay::Real = 0.55,
                   target::Real = 0.45)
    new(logf, A, D, epsilon, decay, target, 0)
  end
end

BIATune(x::Vector; args...) =
   BIATune(x, Nullable{Function}(); args...)

BIATune(x::Vector, logf::Function; args...) =
   BIATune(x, Nullable{Function}(logf); args...)


const BIAVariate = SamplerVariate{BIATune}

function validate(v::BIAVariate)
  n = length(v)
  epsilon = v.tune.epsilon
  0.0 < epsilon < 0.5 ||
    throw(ArgumentError("epsilon is not in (0, 0.5)"))
  length(v.tune.A) == n && all(epsilon .< v.tune.A .< 1 - epsilon) ||
    throw(ArgumentError("A must contain $n probabilities in ($epsilon, $(1 - epsilon))"))
  length(v.tune.D) == n && all(epsilon .< v.tune.D .< 1 - epsilon) ||
    throw(ArgumentError("D must contain $n probabilities in ($epsilon, $(1 - epsilon))"))
  0.5 < v.tune.decay <= 1.0 ||
    throw(ArgumentError("decay is not in (0.5, 1]"))
  0.0 < v.tune.target < 1.0 ||
    throw(ArgumentError("target is not in (0, 1)"))
  validatebinary(v)
end


#################### Sampler Constructor ####################

function BIA(params::ElementOrVector{Symbol}; args...)
  samplerfx = function(model::Model, block::Integer)
    block = SamplingBlock(model, block)
    v = SamplerVariate(block; args...)
    sample!(v, x -> logpdf!(block, x))
    relist(block, v)
  end
  Sampler(params, samplerfx, BIATune())
end


#################### Sampling Functions ####################

sample!(v::BIAVariate) = sample!(v, v.tune.logf)

function sample!(v::BIAVariate, logf::Function)
  tune = v.tune
  tune.iter += 1

  p = length(v)

  ## proposal
  x = v[:]
  added = zeros(Int8, p)
  deleted = zeros(Int8, p)
  q_ratio = 1.0
  for j in 1:p
    if v.value[j] == 0.0
      if rand() < tune.A[j]
        x[j] = 1.0
        added[j] = 1
        q_ratio *= tune.D[j] / tune.A[j]
      end
    else
      if rand() < tune.D[j]
        x[j] = 0.0
        deleted[j] = 1
        q_ratio *= tune.A[j] / tune.D[j]
      end
    end
  end

  ## M-H acceptance probability
  alpha = min(1.0, exp(logf(x) - logf(v.value)) * q_ratio)

  ## Adapt tuning parameters
  for j in 1:p
    ## new A[j]
    C = log((tune.A[j] - tune.epsilon) / (1.0 - tune.A[j] - tune.epsilon)) +
        tune.iter^(-tune.decay) * added[j] * (alpha - tune.target)
    tune.A[j] = (exp(C) * (1.0 - tune.epsilon) + tune.epsilon) / (1.0 + exp(C))

    ## new D[j]
    C = log((tune.D[j] - tune.epsilon) / (1.0 - tune.D[j] - tune.epsilon)) +
        tune.iter^(-tune.decay) * deleted[j] * (alpha - tune.target)
    tune.D[j] = (exp(C) * (1.0 - tune.epsilon) + tune.epsilon) / (1.0 + exp(C))
  end

  ## Accept or reject proposal
  if rand() < alpha
    v[:] = x
  end

  v
end
