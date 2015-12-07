#################### Random Walk Metropolis ####################

#################### Types and Constructors ####################

type RWMTune <: SamplerTune
  delta::Vector{Float64}

  function RWMTune(value::Vector{Float64}=Float64[])
    new(Array(Float64, 0))
  end
end


typealias RWMVariate SamplerVariate{RWMTune}


#################### Sampler Constructor ####################

function RWM{T<:Real}(params::Vector{Symbol}, delta::Vector{T})
  samplerfx = function(model::Model, block::Integer)
    v = SamplerVariate(model, block, true)
    f = x -> logpdf!(model, x, block, true)
    rwm!(v, delta, f)
    relist(model, v, block, true)
  end
  Sampler(params, samplerfx, RWMTune())
end


#################### Sampling Functions ####################

function rwm!{T<:Real}(v::RWMVariate, delta::Vector{T}, logf::Function)
  v.tune.delta = delta

  x = v + delta .* (2.0 * rand(length(v)) - 1.0)

  if rand() < exp(logf(x) - logf(v.value))
    v[:] = x
  end

  v
end
