#################### Metropolis-Adjusted Langevin Algorithm ####################

#################### Types ####################

type MALATune
  scale::Float64
  SigmaF::Cholesky{Float64}
end

type MALAVariate <: VectorVariate
  value::Vector{Float64}
  tune::MALATune

  MALAVariate{T<:Real}(x::AbstractVector{T}, tune::MALATune) = new(x, tune)
end

function MALAVariate{T<:Real}(x::AbstractVector{T}, tune=nothing)
  tune = MALATune(
    NaN,
    Cholesky(Array(Float64, 0, 0), :U)
  )
  MALAVariate(x, tune)
end


#################### Sampler Constructor ####################

function MALA(params::Vector{Symbol}, scale::Real; dtype::Symbol=:forward)
  Sampler(params,
    quote
      tunepar = tune(model, block)
      x = unlist(model, block, true)
      v = MALAVariate(x, tunepar["sampler"])
      fx = x -> malafx!(model, x, block, tunepar["dtype"])
      mala!(v, tunepar["scale"], fx)
      tunepar["sampler"] = v.tune
      relist(model, v, block, true)
    end,
    Dict("scale" => scale, "dtype" => dtype, "sampler" => nothing)
  )
end

function MALA{T<:Real}(params::Vector{Symbol}, scale::Real, Sigma::Matrix{T};
                       dtype::Symbol=:forward)
  Sampler(params,
    quote
      tunepar = tune(model, block)
      x = unlist(model, block, true)
      v = MALAVariate(x, tunepar["sampler"])
      fx = x -> malafx!(model, x, block, tunepar["dtype"])
      mala!(v, tunepar["scale"], tunepar["SigmaF"], fx)
      tunepar["sampler"] = v.tune
      relist(model, v, block, true)
    end,
    Dict("scale" => scale, "SigmaF" => cholfact(Sigma), "dtype" => dtype,
         "sampler" => nothing)
  )
end

function malafx!{T<:Real}(m::Model, x::Vector{T}, block::Integer, dtype::Symbol)
  logf = logpdf!(m, x, block, true)
  grad = isfinite(logf) ?
           gradlogpdf!(m, x, block, true, dtype=dtype) :
           zeros(x)
  logf, grad
end

export malafx!


#################### Sampling Functions ####################

function mala!(v::MALAVariate, scale::Real, fx::Function)
  scale2 = scale / 2.0

  logf0, grad0 = fx(v.value)
  y = v + scale2 * grad0 + sqrt(scale) * randn(length(v))
  logf1, grad1 = fx(y)

  q0 = -0.5 * sumabs2((v - y - scale2 * grad1)) / scale
  q1 = -0.5 * sumabs2((y - v - scale2 * grad0)) / scale

  if rand() < exp((logf1 - q1) - (logf0 - q0))
    v[:] = y
  end
  v.tune.scale = scale

  v
end

function mala!(v::MALAVariate, scale::Real, SigmaF::Cholesky{Float64},
               fx::Function)
  L = sqrt(scale) * SigmaF[:L]
  Linv = inv(L)
  M2 = 0.5 * L * L'

  logf0, grad0 = fx(v.value)
  y = v + M2 * grad0 + L * randn(length(v))
  logf1, grad1 = fx(y)

  q0 = -0.5 * sumabs2(Linv * (v - y - M2 * grad1))
  q1 = -0.5 * sumabs2(Linv * (y - v - M2 * grad0))

  if rand() < exp((logf1 - q1) - (logf0 - q0))
    v[:] = y
  end
  v.tune.scale = scale
  v.tune.SigmaF = SigmaF

  v
end
