#################### Metropolis-Adjusted Langevin Algorithm ####################

#################### Types ####################

type MALATune <: SamplerTune
  logfgrad::Nullable{Function}
  epsilon::Float64
  SigmaL::Union{UniformScaling{Int}, LowerTriangular{Float64}}

  MALATune() = new()

  MALATune(x::Vector, epsilon::Real) = new(Nullable{Function}(), epsilon, I)

  MALATune(x::Vector, epsilon::Real, logfgrad::Function) =
    new(Nullable{Function}(logfgrad), epsilon, I)

  MALATune{T<:Real}(x::Vector, epsilon::Real, Sigma::Matrix{T}) =
    new(Nullable{Function}(), epsilon, cholfact(Sigma)[:L])

  function MALATune{T<:Real}(x::Vector, epsilon::Real, Sigma::Matrix{T},
                             logfgrad::Function)
    new(Nullable{Function}(logfgrad), epsilon, cholfact(Sigma)[:L])
  end
end


const MALAVariate = SamplerVariate{MALATune}

validate(v::MALAVariate) = validate(v, v.tune.SigmaL)

validate(v::MALAVariate, SigmaL::UniformScaling) = v

function validate(v::MALAVariate, SigmaL::LowerTriangular)
  n = length(v)
  size(SigmaL, 1) == n ||
    throw(ArgumentError("Sigma dimension differs from variate length $n"))
  v
end


#################### Sampler Constructor ####################

function MALA(params::ElementOrVector{Symbol}, epsilon::Real; args...)
  MALASampler(params, epsilon; args...)
end

function MALA{T<:Real}(params::ElementOrVector{Symbol}, epsilon::Real,
                       Sigma::Matrix{T}; args...)
  MALASampler(params, epsilon, Sigma; args...)
end

function MALASampler(params, pargs...; dtype::Symbol=:forward)
  samplerfx = function(model::Model, block::Integer)
    block = SamplingBlock(model, block, true)
    v = SamplerVariate(block, pargs...)
    sample!(v, x -> logpdfgrad!(block, x, dtype))
    relist(block, v)
  end
  Sampler(params, samplerfx, MALATune())
end


#################### Sampling Functions ####################

sample!(v::MALAVariate) = sample!(v, v.tune.logfgrad)

function sample!(v::MALAVariate, logfgrad::Function)
  tune = v.tune

  L = sqrt(tune.epsilon) * tune.SigmaL
  Linv = inv(L)
  M2 = 0.5 * L * L'

  logf0, grad0 = logfgrad(v.value)
  y = v + M2 * grad0 + L * randn(length(v))
  logf1, grad1 = logfgrad(y)

  q0 = -0.5 * sumabs2(Linv * (v - y - M2 * grad1))
  q1 = -0.5 * sumabs2(Linv * (y - v - M2 * grad0))

  if rand() < exp((logf1 - q1) - (logf0 - q0))
    v[:] = y
  end

  v
end
