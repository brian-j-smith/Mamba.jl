#################### Random Walk Metropolis ####################

#################### Types and Constructors ####################

type RWMTune <: SamplerTune
  scale::Union{Real, Vector}
  proposal::SymDistributionType

  function RWMTune(value::Vector{Float64}=Float64[])
    new(
      NaN,
      Normal
    )
  end
end


typealias RWMVariate SamplerVariate{RWMTune}


#################### Sampler Constructor ####################

function RWM{T<:Real}(params::Vector{Symbol}, scale::ElementOrVector{T};
                      proposal::SymDistributionType=Normal)
  samplerfx = function(model::Model, block::Integer)
    v = SamplerVariate(model, block, true)
    f = x -> logpdf!(model, x, block, true)
    rwm!(v, scale, f, proposal=proposal)
    relist(model, v, block, true)
  end
  Sampler(params, samplerfx, RWMTune())
end


#################### Sampling Functions ####################

function rwm!{T<:Real}(v::RWMVariate, scale::ElementOrVector{T},
                       logf::Function; proposal::SymDistributionType=Normal)
  v.tune.scale = scale
  v.tune.proposal = proposal

  x = v + scale .* rand(proposal(0.0, 1.0), length(v))

  if rand() < exp(logf(x) - logf(v.value))
    v[:] = x
  end

  v
end
