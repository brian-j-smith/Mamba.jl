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

function HMC(params::Vector{Symbol}, epsilon::Real, L::Int; 
             dtype::Symbol=:forward)
  Sampler(params,
    quote
      tunepar = tune(model, block)
      x = unlist(model, block, true)
      v = HMCVariate(x, tunepar["sampler"])
      fx = x -> hmcfx(model, x, block, tunepar["dtype"])
      hmc!(v, tunepar["epsilon"], tunepar["L"], fx)
      tunepar["sampler"] = v.tune
      relist(model, v, block, true)
    end,
    Dict("epsilon" => epsilon, "L" => L, "dtype" => dtype, "sampler" => nothing)
  )
end

function HMC{T<:Real}(params::Vector{Symbol}, epsilon::Real, L::Int, 
                      Sigma::Matrix{T}; dtype::Symbol=:forward)
  Sampler(params,
    quote
      tunepar = tune(model, block)
      x = unlist(model, block, true)
      v = HMCVariate(x, tunepar["sampler"])
      fx = x -> hmcfx(model, x, block, tunepar["dtype"])
      hmc!(v, tunepar["epsilon"], tunepar["L"], tunepar["SigmaF"], fx)
      tunepar["sampler"] = v.tune
      relist(model, v, block, true)
    end,
    Dict("epsilon" => epsilon, "L" => L, "SigmaF" => cholfact(Sigma), 
         "dtype" => dtype, "sampler" => nothing)
  )
end

function hmcfx{T<:Real}(m::Model, x::Vector{T}, block::Integer, dtype::Symbol)
  logf = logpdf(m, x, block, true)
  grad = isfinite(logf) ?
           gradlogpdf(m, x, block, true, dtype=dtype) :
           zeros(length(x))
  logf, grad
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

function hmc!(v::HMCVariate, epsilon::Float64, L::Int, fx::Function)
  x1 = v[:]
  # Negative log likelihood and gradient
  U0, gradU0 = fx(v)
  U0, gradU0 = -U0, -gradU0
  U1, gradU1 = U0, gradU0

  # momentum variables
  p0 = randn(length(v))
  p1 = p0[:]

  # make a half step for a momentum at the beginning
  p1 = p1 - 0.5 * epsilon * gradU1

  # Alternate full steps for position and momentum
  for i in 1:L
    # Make a full step for the position
    x1 = x1 + epsilon * p1

    U1, gradU1 = fx(x1)
    U1, gradU1 = -U1, -gradU1

    # Make a full step for the momentum, except at the end of trajectory
    if i != L
      p1 = p1 - epsilon * gradU1
    end
  end

  # Make a half step for momentum at the end
  p1 = p1 - 0.5 * epsilon * gradU1

  # Negate momentum at end of trajectory to make the proposal symmetric
  p1 = -p1

  # Evaluate potential and kinetic energies at start and end of trajectory
  Kp0 = 0.5 * sumabs2(p0)
  Kp1 = 0.5 * sumabs2(p1)

  if rand() < exp(U0 - U1 + Kp0 - Kp1)
    v[:] = x1
  end
  v.tune.epsilon = epsilon
  v.tune.L = L

  v
end

function hmc!(v::HMCVariate, epsilon::Float64, L::Int, 
              SigmaF::Cholesky{Float64}, fx::Function)
  S = SigmaF[:L]
  Sinv = inv(L)

  x1 = v[:]
  # Negative log likelihood and gradient
  U0, gradU0 = fx(v)
  U0, gradU0 = -U0, -gradU0
  U1, gradU1 = U0, gradU0

  # momentum variables
  p0 = S * randn(length(v))
  p1 = p0[:]

  # make a half step for a momentum at the beginning
  p1 = p1 - 0.5 * epsilon * gradU1

  # Alternate full steps for position and momentum
  for i in 1:L
    # Make a full step for the position
    x1 = x1 + epsilon * p1

    U1, gradU1 = fx(x1)
    U1, gradU1 = -U1, -gradU1

    # Make a full step for the momentum, except at the end of trajectory
    if i != L
      p1 = p1 - epsilon * gradU1
    end
  end

  # Make a half step for momentum at the end
  p1 = p1 - 0.5 * epsilon * gradU1

  # Negate momentum at end of trajectory to make the proposal symmetric
  p1 = -p1

  # Evaluate potential and kinetic energies at start and end of trajectory
  Kp0 = 0.5 * sumabs2(Sinv * p0)
  Kp1 = 0.5 * sumabs2(Sinv * p1)

  if rand() < exp(U0 - U1 + Kp0 - Kp1)
    v[:] = x1
  end
  v.tune.epsilon = epsilon
  v.tune.L = L
  v.tune.SigmaF = SigmaF

  v
end
