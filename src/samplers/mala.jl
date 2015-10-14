#################### Metropolis Adjusted Langevin Algorithm ####################

#################### Types ####################

type MALATune
  U::Cholesky{Float64}
  M::Matrix{Float64}
  scale::Float64
end

type MALAVariate <: VectorVariate
  value::Vector{Float64}
  tune::MALATune

  MALAVariate(x::Vector{Float64}, tune::MALATune) = new(x, tune)
end

function MALAVariate(x::Vector{Float64}, tune=nothing)
  tune = MALATune(
    Cholesky(Array(Float64, 0, 0), :U),
    Array(Float64, 0, 0),
    NaN
  )
  MALAVariate(x, tune)
end


#################### Sampler Constructor ####################

function MALA(params::Vector{Symbol}, scale::Float64, M::Matrix{Float64}, dtype::Symbol=:forward)
  Sampler(params,
    quote
      tunepar = tune(model, block)
      x = unlist(model, block, true)
      v = MALAVariate(x, tunepar["sampler"])
      fx = x -> malafx!(model, x, block, tunepar["dtype"])

      mala!(v, tunepar["scale"], tunepar["U"], fx)
      tunepar["sampler"] = v.tune
      relist(model, v.value, block, true)
    end,
    Dict("scale" => scale, "U" => cholfact(M), "dtype" => dtype, "sampler" => nothing)
  )
end

function malafx{T<:Real}(m::Model, x::Vector{T}, block::Integer, dtype::Symbol)
  logf = logpdf(m, x, block, true)
  grad = isfinite(logf) ?
           gradlogpdf(m, x, block, true, dtype=dtype) :
           zeros(length(x))
  logf, grad
end

function malafx!{T<:Real}(m::Model, x::Vector{T}, block::Integer, dtype::Symbol)
  logf = logpdf!(m, x, block, true)
  grad = isfinite(logf) ?
           gradlogpdf!(m, x, block, true, dtype=dtype) :
           zeros(length(x))
  logf, grad
end

export malafx, malafx!

#################### Sampling Functions ####################

function mala!(v::MALAVariate, scale::Float64, U::Cholesky{Float64}, fx::Function)
  v.tune.scale = scale
  v.tune.U = U
  M = U[:L]*U[:L]'

  logf0, grad0 = fx(v.value)
  y = v + 0.5 * scale^2 * M * grad0 + scale * U[:L] * randn(length(v))

  logf1, grad1 = fx(y)

  q0 = logpdf(MvNormal(v + 0.5 * scale^2 * M * grad0, scale^2 * M), y)
  q1 = logpdf(MvNormal(y + 0.5 * scale^2 * M * grad1, scale^2 * M), v)

  if rand() < exp( (logf1 + q1) - (logf0 + q0))
    v[:] = y
  end
  v
end
