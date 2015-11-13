############ Binary Gibbs Sampler ##############

#################### Types ####################

type BGSTune
  probs::Vector{Float64}
end

type BGSVariate <: VectorVariate
  value::Vector{Float64}
  tune::BGSTune

  function BGSVariate(x::Vector{Float64}, tune::BGSTune)
    all(insupport(Bernoulli, x)) ||
      throw(ArgumentError("x is not a binary vector"))
    new(x, tune)
  end
end

function BGSVariate(x::Vector{Float64}, tune=nothing)
  tune = BGSTune(
    fill(NaN, length(x))
  )
  BGSVariate(x, tune)
end


#################### Sampler Constructor ####################

function BGS(params::Vector{Symbol})
  Sampler(params,
    quote
      tunepar = tune(model, block)
      x = unlist(model, block)
      v = BGSVariate(x, tunepar["sampler"])
      f = x -> logpdf!(model, x, block)
      bgs!(v, f)
      tunepar["sampler"] = v.tune
      relist(model, v, block)
    end,
    Dict{AbstractString,Any}("sampler" => nothing)
  )
end


#################### Sampling Functions ####################

function bgs!(v::BGSVariate, logf::Function)
  tune = v.tune

  for i in 1:length(v)
    v0 = v[:]
    v1 = v[:]

    v0[i] = 0.0
    v1[i] = 1.0
    p = invlogit(logf(v1) - logf(v0))
    if p < 0.0 || p > 1.0
      p = 0.5
    end
    v[i] = rand() < p ? 1.0 : 0.0
    tune.probs[i] = p
  end
  v
end
