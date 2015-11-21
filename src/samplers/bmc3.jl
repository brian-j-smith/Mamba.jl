#################### Binary MCMC Model Composition ####################

#################### Types and Constructors ####################

type BMC3Tune
  indexset::Vector{Vector{Int}}
end

type BMC3Variate <: VectorVariate
  value::Vector{Float64}
  tune::BMC3Tune

  function BMC3Variate{T<:Real}(x::AbstractVector{T}, tune::BMC3Tune)
    all(insupport(Bernoulli, x)) ||
      throw(ArgumentError("x is not a binary vector"))
    new(x, tune)
  end
end

function BMC3Variate{T<:Real}(x::AbstractVector{T}, tune=nothing)
  tune = BMC3Tune(
    Vector{Vector{Int}}[]
  )
  BMC3Variate(x, tune)
end


#################### Sampler Constructor ####################

function BMC3(params::Vector{Symbol}; k::Integer=1)
  Sampler(params, (model::Model, block::Integer) ->
    begin
      tunepar = tune(model, block)
      x = unlist(model, block)
      v = BMC3Variate(x, tunepar["sampler"])
      f = x -> logpdf!(model, x, block)
      bmc3!(v, f, k=tunepar["k"])
      tunepar["sampler"] = v.tune
      relist(model, v, block)
    end,
    Dict("k" => k, "sampler" => nothing)
  )
end

function BMC3(params::Vector{Symbol}, indexset::Vector{Vector{Int}})
  Sampler(params, (model::Model, block::Integer) ->
    begin
      tunepar = tune(model, block)
      x = unlist(model, block)
      v = BMC3Variate(x, tunepar["sampler"])
      f = x -> logpdf!(model, x, block)
      bmc3!(v, tunepar["indexset"], f)
      tunepar["sampler"] = v.tune
      relist(model, v, block)
    end,
    Dict("indexset" => indexset, "sampler" => nothing)
  )
end


#################### Sampling Functions ####################

function bmc3!(v::BMC3Variate, logf::Function; k::Integer=1)
  d = length(v)
  k <= d || throw(ArgumentError("k is greater than length(v)"))
  bmc3_sub!(v, randperm(d)[1:k], logf)
end

function bmc3!(v::BMC3Variate, indexset::Vector{Vector{Int}}, logf::Function)
  v.tune.indexset = indexset
  bmc3_sub!(v, indexset[rand(1:length(indexset))], logf)
end

function bmc3_sub!(v::BMC3Variate, idx::Vector{Int}, logf::Function)
  x = v[:]
  x[idx] = 1.0 - v[idx]
  if rand() < exp(logf(x) - logf(v.value))
    v[:] = x
  end
  v
end
