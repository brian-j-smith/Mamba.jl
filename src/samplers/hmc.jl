#################### Hamiltonian Monte Carlo ####################

#################### Types and Constructors ####################

type HMCTune <: SamplerTune
  logfgrad::Nullable{Function}
  epsilon::Float64
  L::Int
  SigmaL::Union{UniformScaling{Int}, LowerTriangular{Float64}}

  HMCTune() = new()

  HMCTune(x::Vector, epsilon::Real, L::Integer) =
    new(Nullable{Function}(), epsilon, L, I)

  HMCTune(x::Vector, epsilon::Real, L::Integer, logfgrad::Function) =
    new(Nullable{Function}(logfgrad), epsilon, L, I)

  function HMCTune{T<:Real}(x::Vector, epsilon::Real, L::Integer,
                            Sigma::Matrix{T})
    new(Nullable{Function}(), epsilon, L, cholfact(Sigma)[:L])
  end

  function HMCTune{T<:Real}(x::Vector, epsilon::Real, L::Integer,
                            Sigma::Matrix{T}, logfgrad::Function)
    new(Nullable{Function}(logfgrad), epsilon, L, cholfact(Sigma)[:L])
  end
end


typealias HMCVariate SamplerVariate{HMCTune}

validate(v::HMCVariate) = validate(v, v.tune.SigmaL)

validate(v::HMCVariate, SigmaL::UniformScaling) = v

function validate(v::HMCVariate, SigmaL::LowerTriangular)
  n = length(v)
  size(SigmaL, 1) == n ||
    throw(ArgumentError("Sigma dimension differs from variate length $n"))
  v
end


#################### Sampler Constructor ####################

function HMC(params::ElementOrVector{Symbol}, epsilon::Real, L::Integer;
             dtype::Symbol=:forward)
  HMCSampler(params, dtype, epsilon, L)
end

function HMC{T<:Real}(params::ElementOrVector{Symbol}, epsilon::Real,
                      L::Integer, Sigma::Matrix{T}; dtype::Symbol=:forward)
  HMCSampler(params, dtype, epsilon, L, Sigma)
end

function HMCSampler(params, dtype, pargs...)
  samplerfx = function(model::Model, block::Integer)
    block = SamplingBlock(model, block, true)
    v = SamplerVariate(block, pargs...)
    sample!(v, x -> logpdfgrad!(block, x, dtype))
    relist(block, v)
  end
  Sampler(params, samplerfx, HMCTune())
end


#################### Sampling Functions ####################

sample!(v::HMCVariate) = sample!(v, v.tune.logfgrad)

function sample!(v::HMCVariate, logfgrad::Function)
  tune = v.tune

  x1 = v[:]
  logf0, grad0 = logf1, grad1 = logfgrad(x1)

  ## Momentum variables
  p0 = p1 = tune.SigmaL * randn(length(v))

  ## Make a half step for a momentum at the beginning
  p1 += 0.5 * tune.epsilon * grad0

  ## Alternate full steps for position and momentum
  for i in 1:tune.L
    ## Make a full step for the position
    x1 += tune.epsilon * p1

    logf1, grad1 = logfgrad(x1)

    ## Make a full step for the momentum
    p1 += tune.epsilon * grad1
  end

  ## Make a half step for momentum at the end
  p1 -= 0.5 * tune.epsilon * grad1

  ## Negate momentum at end of trajectory to make the proposal symmetric
  p1 *= -1.0

  ## Evaluate potential and kinetic energies at start and end of trajectory
  SigmaLinv = inv(tune.SigmaL)
  Kp0 = 0.5 * sumabs2(SigmaLinv * p0)
  Kp1 = 0.5 * sumabs2(SigmaLinv * p1)

  if rand() < exp((logf1 - Kp1) - (logf0 - Kp0))
    v[:] = x1
  end

  v
end
