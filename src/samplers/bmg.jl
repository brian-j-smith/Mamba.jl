############ Binary Metropolised Gibbs Sampler ##############

#################### Types ####################

type BMGTune
  probs::Vector{Float64}
end

type BMGVariate <: VectorVariate
  value::Vector{Float64}
  tune::BMGTune

  function BMGVariate{T<:Real}(x::AbstractVector{T}, tune::BMGTune)
    all(insupport(Bernoulli, x)) ||
      throw(ArgumentError("x is not a binary vector"))
    new(x, tune)
  end
end

function BMGVariate{T<:Real}(x::AbstractVector{T}, tune=nothing)
  tune = BMGTune(
    fill(NaN, length(x))
  )
  BMGVariate(x, tune)
end


#################### Sampler Constructor ####################

function BMG(params::Vector{Symbol})
  Sampler(params,
    quote
      tunepar = tune(model, block)
      x = unlist(model, block)
      v = BMGVariate(x, tunepar["sampler"])
      f = x -> logpdf!(model, x, block)
      bmg!(v, f)
      tunepar["sampler"] = v.tune
      relist(model, v, block)
    end,
    Dict{AbstractString,Any}("sampler" => nothing)
  )
end


#################### Sampling Functions ####################

function bmg!(v::BMGVariate, logf::Function)
  y = v[:]
  tune = v.tune

  v0 = v[:]
  v1 = v[:]
  for i in 1:length(v)
    v0_i = v0[i]
    v1_i = v1[i]

    v0[i] = 0.0
    v1[i] = 1.0
    p = invlogit(logf(v1) - logf(v0))
    if p < 0.0 || p > 1.0
      p = 0.5
    end
    y[i] = rand() < p ? 1.0 : 0.0
    tune.probs[i] = p

    v0[i] = v0_i
    v1[i] = v1_i
  end
  if rand() < exp(logf(y) - logf(v.value))
    v[:] = y
  end

  v
end
