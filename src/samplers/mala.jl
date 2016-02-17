#################### Metropolis-Adjusted Langevin Algorithm ####################

#################### Types ####################

type MALATune <: SamplerTune
  logfgrad::Nullable{Function}
  scale::Float64
  SigmaL::Union{UniformScaling{Int}, LowerTriangular{Float64}}

  MALATune() = new()

  MALATune(x::Vector, scale::Real) = new(Nullable{Function}(), scale, I)

  MALATune(x::Vector, scale::Real, logfgrad::Function) =
    new(Nullable{Function}(logfgrad), scale, I)

  MALATune{T<:Real}(x::Vector, scale::Real, Sigma::Matrix{T}) =
    new(Nullable{Function}(), scale, cholfact(Sigma)[:L])

  function MALATune{T<:Real}(x::Vector, scale::Real, Sigma::Matrix{T},
                             logfgrad::Function)
    new(Nullable{Function}(logfgrad), scale, cholfact(Sigma)[:L])
  end
end


typealias MALAVariate SamplerVariate{MALATune}

validate(v::MALAVariate) = validate(v, v.tune.SigmaL)

validate(v::MALAVariate, SigmaL::UniformScaling) = v

function validate(v::MALAVariate, SigmaL::LowerTriangular)
  n = length(v)
  size(SigmaL, 1) == n ||
    throw(ArgumentError("Sigma dimension differs from variate length $n"))
  v
end


#################### Sampler Constructor ####################

function MALA(params::ElementOrVector{Symbol}, scale::Real;
              dtype::Symbol=:forward)
  MALASampler(params, dtype, scale)
end

function MALA{T<:Real}(params::ElementOrVector{Symbol}, scale::Real,
                       Sigma::Matrix{T}; dtype::Symbol=:forward)
  MALASampler(params, dtype, scale, Sigma)
end

function MALASampler(params, dtype, pargs...)
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

  L = sqrt(tune.scale) * tune.SigmaL
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
