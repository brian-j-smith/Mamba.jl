#################### Hamiltonian Monte Carlo ####################

#################### Types and Constructors ####################

type HMCTune
  epsilon::Float64
  L::Int
  SigmaF::Cholesky{Float64}
end

type HMCVariate <: VectorVariate
  value::Vector{Float64}
  tune::HMCTune

  HMCVariate(x::Vector{Float64}, tune::HMCTune) = new(x, tune)
end

function HMCVariate(x::Vector{Float64}, tune=nothing)
  tune = HMCTune(
    NaN,
    0,
    Cholesky(Array(Float64, 0, 0), :U)
  )
  HMCVariate(x, tune)
end


#################### Sampler Constructor ####################

function HMC(params::Vector{Symbol}, epsilon::Real, L::Integer;
             dtype::Symbol=:forward)
  Sampler(params,
    quote
      tunepar = tune(model, block)
      x = unlist(model, block, true)
      v = HMCVariate(x, tunepar["sampler"])
      fx = x -> hmcfx!(model, x, block, tunepar["dtype"])
      hmc!(v, tunepar["epsilon"], tunepar["L"], fx)
      tunepar["sampler"] = v.tune
      relist(model, v, block, true)
    end,
    Dict("epsilon" => epsilon, "L" => L, "dtype" => dtype, "sampler" => nothing)
  )
end

function HMC{T<:Real}(params::Vector{Symbol}, epsilon::Real, L::Integer,
                      Sigma::Matrix{T}; dtype::Symbol=:forward)
  Sampler(params,
    quote
      tunepar = tune(model, block)
      x = unlist(model, block, true)
      v = HMCVariate(x, tunepar["sampler"])
      fx = x -> hmcfx!(model, x, block, tunepar["dtype"])
      hmc!(v, tunepar["epsilon"], tunepar["L"], tunepar["SigmaF"], fx)
      tunepar["sampler"] = v.tune
      relist(model, v, block, true)
    end,
    Dict("epsilon" => epsilon, "L" => L, "SigmaF" => cholfact(Sigma),
         "dtype" => dtype, "sampler" => nothing)
  )
end

function hmcfx!{T<:Real}(m::Model, x::Vector{T}, block::Integer, dtype::Symbol)
  logf = logpdf!(m, x, block, true)
  grad = isfinite(logf) ?
           gradlogpdf!(m, x, block, true, dtype=dtype) :
           zeros(length(x))
  logf, grad
end

export hmcfx!


#################### Sampling Functions ####################

function hmc!(v::HMCVariate, epsilon::Real, L::Integer, fx::Function)
  x1 = v[:]
  logf0, grad0 = logf1, grad1 = fx(x1)

  ## Momentum variables
  p0 = p1 = randn(length(v))

  ## Make a half step for a momentum at the beginning
  p1 += 0.5 * epsilon * grad0

  ## Alternate full steps for position and momentum
  for i in 1:L
    ## Make a full step for the position
    x1 += epsilon * p1

    logf1, grad1 = fx(x1)

    ## Make a full step for the momentum
    p1 += epsilon * grad1
  end

  ## Make a half step for momentum at the end
  p1 -= 0.5 * epsilon * grad1

  ## Negate momentum at end of trajectory to make the proposal symmetric
  p1 *= -1.0

  ## Evaluate potential and kinetic energies at start and end of trajectory
  Kp0 = 0.5 * sumabs2(p0)
  Kp1 = 0.5 * sumabs2(p1)

  if rand() < exp((logf1 - Kp1) - (logf0 - Kp0))
    v[:] = x1
  end
  v.tune.epsilon = epsilon
  v.tune.L = L

  v
end

function hmc!(v::HMCVariate, epsilon::Real, L::Integer,
              SigmaF::Cholesky{Float64}, fx::Function)
  S = SigmaF[:L]
  Sinv = inv(S)

  x1 = v[:]
  logf0, grad0 = logf1, grad1 = fx(x1)

  ## Momentum variables
  p0 = p1 = S * randn(length(v))

  ## Make a half step for a momentum at the beginning
  p1 += 0.5 * epsilon * grad0

  ## Alternate full steps for position and momentum
  for i in 1:L
    ## Make a full step for the position
    x1 += epsilon * p1

    logf1, grad1 = fx(x1)

    ## Make a full step for the momentum
    p1 += epsilon * grad1
  end

  ## Make a half step for momentum at the end
  p1 -= 0.5 * epsilon * grad1

  ## Negate momentum at end of trajectory to make the proposal symmetric
  p1 *= -1.0

  ## Evaluate potential and kinetic energies at start and end of trajectory
  Kp0 = 0.5 * sumabs2(Sinv * p0)
  Kp1 = 0.5 * sumabs2(Sinv * p1)

  if rand() < exp((logf1 - Kp1) - (logf0 - Kp0))
    v[:] = x1
  end
  v.tune.epsilon = epsilon
  v.tune.L = L
  v.tune.SigmaF = SigmaF

  v
end
