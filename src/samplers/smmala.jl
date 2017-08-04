########## Simplified Manifold Metropolis-Adjusted Langevin Algorithm ##########

#################### Types ####################

type SMMALATune <: SamplerTune
  logfgradhess::Nullable{Function}
  epsilon::Float64
  SigmaL::Union{UniformScaling{Int}, LowerTriangular{Float64}}

  SMMALATune() = new()

  SMMALATune(x::Vector, epsilon::Real) = new(Nullable{Function}(), epsilon, I)

  SMMALATune(x::Vector, epsilon::Real, logfgradhess::Function) =
    new(Nullable{Function}(logfgradhess), epsilon, I)

  SMMALATune{T<:Real}(x::Vector, epsilon::Real, Sigma::Matrix{T}) =
    new(Nullable{Function}(), epsilon, cholfact(Sigma)[:L])

  function SMMALATune{T<:Real}(x::Vector, epsilon::Real, Sigma::Matrix{T},
                             logfgradhess::Function)
    new(Nullable{Function}(logfgradhess), epsilon, cholfact(Sigma)[:L])
  end
end


const SMMALAVariate = SamplerVariate{SMMALATune}

validate(v::SMMALAVariate) = validate(v, v.tune.SigmaL)

validate(v::SMMALAVariate, SigmaL::UniformScaling) = v

function validate(v::SMMALAVariate, SigmaL::LowerTriangular)
  n = length(v)
  size(SigmaL, 1) == n ||
    throw(ArgumentError("Sigma dimension differs from variate length $n"))
  v
end


#################### Sampler Constructor ####################

function SMMALA(params::ElementOrVector{Symbol}, epsilon::Real; args...)
  SMMALASampler(params, epsilon; args...)
end

function SMMALA{T<:Real}(params::ElementOrVector{Symbol}, epsilon::Real,
                       Sigma::Matrix{T}; args...)
  SMMALASampler(params, epsilon, Sigma; args...)
end

function SMMALASampler(params, pargs...)
  samplerfx = function(model::Model, block::Integer)
    block = SamplingBlock(model, block, true)
    v = SamplerVariate(block, pargs...)
    sample!(v, x -> logpdfgradhess!(block, x))
    relist(block, v)
  end
  Sampler(params, samplerfx, SMMALATune())
end


#################### Sampling Functions ####################

sample!(v::SMMALAVariate) = sample!(v, v.tune.logfgradhess)

function sample!(v::SMMALAVariate, logfgradhess::Function)
  tune = v.tune

  L = sqrt(tune.epsilon) * tune.SigmaL
  Linv = inv(L)
  M2 = 0.5 * L * L'

  logf0, grad0, hess0 = logfgradhess(v.value)
  logf1, grad1, hess1 = logfgradhess(y)

  sqrthess0 = cholfact(hess0)[:L]
  sqrthess1 = cholfact(hess1)[:L]

  invhess0 = inv(hess0)
  invhess1 = inv(hess1)

  sqrtinvhess0 = cholfact(Hermitian(invhess0))[:L]
  sqrtinvhess1 = cholfact(Hermitian(invhess1))[:L]

  y = v + M2 * invhess0 * grad0 + sqrtinvhess0 * L * randn(length(v))

  q0 = -0.5 * sum(abs2, (Linv * sqrthess1) * (v - y - M2 * invhess1 * grad1))
  q1 = -0.5 * sum(abs2, (Linv * sqrthess0) * (y - v - M2 * invhess0 * grad0))

  if rand() < exp((logf1 - q1) - (logf0 - q0))
    v[:] = y
  end

  v
end
