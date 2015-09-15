#################### Binary Deterministic Sampler ####################

#################### Types ####################

type BDSTune
  indexset::Vector{Vector{Int}}
end

type BDSVariate <: VectorVariate
  value::Vector{Float64}
  tune::BDSTune

  BDSVariate(x::Vector{Float64}, tune::BDSTune) = new(x, tune)
end

function BDSVariate(x::Vector{Float64}, tune=nothing)
  tune = BDSTune(
    Vector{Vector{Int}}[]
  )
  BDSVariate(x, tune)
end


#################### Sampler Constructor ####################

function BDS(params::Vector{Symbol}, d::Integer, k::Integer=1)
  d >= k > 0 || throw(ArgumentError("values must be d >= k > 0"))
  indexset = collect(combinations(1:d, k))
  BDS(params, indexset)
end

function BDS(params::Vector{Symbol}, indexset::Vector{Vector{Int}})
  Sampler(params,
    quote
      x = unlist(model, block, false)
      tunepar = tune(model, block)
      f = y -> logpdf!(model, y, block, false)
      v = BDSVariate(x)
      bds!(v, tunepar["indexset"], f)
      relist(model, v.value, block, false)
    end,
    Dict("indexset" => indexset)
  )
end


#################### Sampling Functions ####################

function bds!(v::BDSVariate, indexset::Vector{Vector{Int}}, logf::Function)
  x = v[:]

  # sample indices and flip values
  idx = indexset[rand(1:length(indexset))]
  x[idx] = 1.0 - v[idx]

  if rand() < exp(logf(x) - logf(v.value))
    v[:] = x
  end
  v.tune.indexset = indexset

  v
end
