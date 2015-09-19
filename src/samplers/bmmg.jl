#################### Binary Modified Metropolised Gibbs ####################

#################### Types ####################

type BMMGTune
  indexset::Vector{Vector{Int}}
end

type BMMGVariate <: VectorVariate
  value::Vector{Float64}
  tune::BMMGTune

  BMMGVariate(x::Vector{Float64}, tune::BMMGTune) = new(x, tune)
end

function BMMGVariate(x::Vector{Float64}, tune=nothing)
  tune = BMMGTune(
    Vector{Vector{Int}}[]
  )
  BMMGVariate(x, tune)
end


#################### Sampler Constructor ####################

function BMMG(params::Vector{Symbol}, d::Integer, k::Integer=1)
  d >= k > 0 || throw(ArgumentError("values must be d >= k > 0"))
  indexset = collect(combinations(1:d, k))
  BMMG(params, indexset)
end

function BMMG(params::Vector{Symbol}, indexset::Vector{Vector{Int}})
  Sampler(params,
    quote
      x = unlist(model, block, false)
      tunepar = tune(model, block)
      f = y -> logpdf!(model, y, block, false)
      v = BMMGVariate(x)
      bmmg!(v, tunepar["indexset"], f)
      relist(model, v.value, block, false)
    end,
    Dict("indexset" => indexset)
  )
end


#################### Sampling Functions ####################

function bmmg!(v::BMMGVariate, indexset::Vector{Vector{Int}}, logf::Function)
  all(insupport(Bernoulli, v)) ||
    throw(ArgumentError("must supply a binary vector"))

  x = v[:]

  # sample indices and change values
  idx = indexset[rand(1:length(indexset))]
  x[idx] = 1.0 - v[idx]

  if rand() < exp(logf(x) - logf(v.value))
    v[:] = x
  end
  v.tune.indexset = indexset

  v
end
