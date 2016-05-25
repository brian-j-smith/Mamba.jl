#################### Binary MCMC Model Composition ####################

#################### Types and Constructors ####################

type BMC3Tune <: SamplerTune
  logf::Nullable{Function}
  k::Int
  indexset::Vector{Vector{Int}}

  n::Int
  m::Vector{Float64}
  v::Vector{Float64}

  epsilon::Float64

  BMC3Tune() = new()

  BMC3Tune(x::Vector, logf::Nullable{Function}; k::Integer=1, 
           indexset::Vector{Vector{Int}} = Vector{Vector{Int}}(),
           epsilon::Float64 = 1/length(x)) =
    new(logf, k, indexset, 0, zeros(length(x)), zeros(length(x)), epsilon)
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

function BMC3(params::ElementOrVector{Symbol}; adapt::Symbol = :none, args...)

  adapt in [:var, :freq, :none] || 
    throw(ArgumentError("adapt must be one of :var, :freq, or :none"))
    
  samplerfx = function(model::Model, block::Integer)
    block = SamplingBlock(model, block)
    v = SamplerVariate(block; args...)
    sample!(v, x -> logpdf!(block, x), adapt=adapt)
    relist(block, v)
  end
  Sampler(params, samplerfx, BMC3Tune())
end


#################### Sampling Functions ####################

sample!(v::BMC3Variate; args...) = sample!(v, v.tune.logf; args...)

function sample!(v::BMC3Variate, logf::Function; adapt::Symbol = :none)
  x = v[:]

  ## online calculation of mean and variance
  v.tune.n += 1
  d = x - v.tune.m
  v.tune.m += d/v.tune.n
  v.tune.v += d .* (x - v.tune.m)

  idx = Int64[]
  if adapt == :var
    idx = wsample(1:length(v), 
            (1 - v.tune.epsilon) * (v.tune.v / (v.tune.n - 1)) + v.tune.epsilon)
  elseif adapt == :freq
    idx = wsample(1:length(v), (1 - v.tune.epsilon) * v.tune.m + v.tune.epsilon)
  else
    if length(v.tune.indexset) > 0
      idx = rand(v.tune.indexset)
    else
      idx = randperm(length(v))[1:v.tune.k]
    end
  end
  x[idx] = 1.0 - v[idx]
  if rand() < exp(logf(x) - logf(v.value))
    v[:] = x
  end
  v
end
