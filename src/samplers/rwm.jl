#################### Random Walk Metropolis ####################

#################### Types and Constructors ####################

type RWMTune <: SamplerTune
  logf::Nullable{Function}
  scale::Union{Float64, Vector{Float64}}
  proposal::SymDistributionType

  RWMTune() = new()

  function RWMTune(x::Vector, scale::Real, logf::Nullable{Function};
                   proposal::SymDistributionType=Normal)
    new(logf, Float64(scale), proposal)
  end

  function RWMTune{T<:Real}(x::Vector, scale::Vector{T},
                            logf::Nullable{Function};
                            proposal::SymDistributionType=Normal)
    new(logf, convert(Vector{Float64}, scale), proposal)
  end
end

RWMTune{T<:Real}(x::Vector, scale::ElementOrVector{T}; args...) =
  RWMTune(x, scale, Nullable{Function}(); args...)

function RWMTune{T<:Real}(x::Vector, scale::ElementOrVector{T}, logf::Function;
                          args...)
  RWMTune(x, scale, Nullable{Function}(logf); args...)
end


typealias RWMVariate SamplerVariate{RWMTune}

validate(v::RWMVariate) = validate(v, v.tune.scale)

validate(v::RWMVariate, scale::Float64) = v

function validate(v::RWMVariate, scale::Vector)
  n = length(v)
  length(scale) == n ||
    throw(ArgumentError("length(scale) differs from variate length $n"))
  v
end


#################### Sampler Constructor ####################

function RWM{T<:Real}(params::ElementOrVector{Symbol},
                      scale::ElementOrVector{T}; args...)
  samplerfx = function(model::Model, block::Integer)
    block = SamplingBlock(model, block, true)
    v = SamplerVariate(block, scale; args...)
    sample!(v, x -> logpdf!(block, x))
    relist(block, v)
  end
  Sampler(params, samplerfx, RWMTune())
end


#################### Sampling Functions ####################

sample!(v::RWMVariate) = sample!(v, v.tune.logf)

function sample!(v::RWMVariate, logf::Function)
  x = v + v.tune.scale .* rand(v.tune.proposal(0.0, 1.0), length(v))
  if rand() < exp(logf(x) - logf(v.value))
    v[:] = x
  end
  v
end


#################### Legacy Sampler Code ####################

RWMTune(x) = RWMTune(x, NaN)


function rwm!{T<:Real}(v::RWMVariate, scale::ElementOrVector{T},
                       logf::Function; proposal::SymDistributionType=Normal)
  v.tune.scale = scale
  v.tune.proposal = proposal
  sample!(v, logf)
end
