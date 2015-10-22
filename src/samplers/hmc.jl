#################### Hamiltonian Monte Carlo ####################

#################### Types and Constructors ####################

type HMCTune
  M::Matrix{Float64}
  epsilon::Float64
  L::Int
end

type HMCVariate <: VectorVariate
  value::Vector{Float64}
  tune::HMCTune

  HMCVariate(x::Vector{Float64}, tune::HMCTune) = new(x, tune)
end

function HMCVariate(x::Vector{Float64}, tune=nothing)
  tune = HMCTune(
    Array(Float64, 0, 0),
    NaN,
    0
  )
  HMCVariate(x, tune)
end


#################### Sampler Constructor ####################

function HMC{T<:Real}(params::Vector{Symbol}, M::Matrix{T}, epsilon::Float64, L::Int; dtype::Symbol=:forward)
  Sampler(params,
    quote
      tunepar = tune(model, block)
      x = unlist(model, block, true)
      v = HMCVariate(x, tunepar["sampler"])
      f = x -> hmcfx(model, x, block, tunepar["dtype"])
      hmc!(v, tunepar["M"], tunepar["epsilon"], tunepar["L"], f)
      tunepar["sampler"] = v.tune
      relist(model, v, block, true)
    end,
    Dict("M" => M, "epsilon" => epsilon, "L" => L, "dtype" => dtype, "sampler" => nothing)
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

export hmcfx, hmcfx!

#################### Sampling Functions ####################

function hmc!(v::HMCVariate, M::Matrix{Float64}, epsilon::Float64, L::Int, fx::Function)
  tune = v.tune
  tune.M = M
  tune.epsilon = epsilon
  tune.L = L
  
  x1 = v[:]
  logf0, grad0 = fx(v)
  logf1, grad1 = logf0, grad0

  # momentum variables
  p0 = cholfact(M)[:L]*randn(length(v))
  p1 = p0[:]

  # make a half step for a momentum at the beginning
  p1 = p1 + epsilon * grad1 / 2 

  # Alternate full steps for position and momentum
  for i in 1:L
    # Make a full step for the position
    x1 = x1 + epsilon * p1

    logf1, grad1 = fx(x1)

    # Make a full step for the momentum, except at the end of trajectory
    if i != L
      p1 = p1 + epsilon * grad1
    end
  end

  # Make a half step for momentum at the end
  p1 = p1 + epsilon * grad1 / 2

  # Negate momentum at end of trajectory to make the proposal symmetric
  p1 = -p1

  Kp0 = (p0' * inv(M) * p0)[1] /2
  Kp1 = (p1' * inv(M) * p1)[1] /2

  if rand() < exp(logf1 - logf0 + Kp0 - Kp1)
    v[:] = x1
  end
  v
end
