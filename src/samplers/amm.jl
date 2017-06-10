#################### Adaptive Mixture Metropolis ####################

#################### Types and Constructors ####################

type AMMTune <: SamplerTune
  logf::Nullable{Function}
  adapt::Bool
  beta::Float64
  m::Int
  Mv::Vector{Float64}
  Mvv::Matrix{Float64}
  scale::Float64
  SigmaL::LowerTriangular{Float64}
  SigmaLm::Matrix{Float64}

  AMMTune() = new()

  function AMMTune{T<:Real}(x::Vector, Sigma::Matrix{T},
                            logf::Nullable{Function}; beta::Real=0.05,
                            scale::Real=2.38)
    new(logf, false, beta, 0, Vector{Float64}(), Matrix{Float64}(), scale,
        cholfact(Sigma)[:L], Matrix{Float64}())
  end
end

AMMTune{T<:Real}(x::Vector, Sigma::Matrix{T}; args...) =
  AMMTune(x, Sigma, Nullable{Function}(); args...)

AMMTune{T<:Real}(x::Vector, Sigma::Matrix{T}, logf::Function; args...) =
  AMMTune(x, Sigma, Nullable{Function}(logf); args...)


const AMMVariate = SamplerVariate{AMMTune}

function validate(v::AMMVariate)
  n = length(v)
  size(v.tune.SigmaL, 1) == n ||
    throw(ArgumentError("Sigma dimension differs from variate length $n"))
  v
end


#################### Sampler Constructor ####################

function AMM{T<:Real}(params::ElementOrVector{Symbol}, Sigma::Matrix{T};
                      adapt::Symbol=:all, args...)
  adapt in [:all, :burnin, :none] ||
    throw(ArgumentError("adapt must be one of :all, :burnin, or :none"))

  samplerfx = function(model::Model, block::Integer)
    block = SamplingBlock(model, block, true)
    v = SamplerVariate(block, Sigma; args...)
    isadapt = adapt == :burnin ? model.iter <= model.burnin :
              adapt == :all ? true : false
    sample!(v, x -> logpdf!(block, x), adapt=isadapt)
    relist(block, v)
  end
  Sampler(params, samplerfx, AMMTune())
end


#################### Sampling Functions ####################

sample!(v::AMMVariate; args...) = sample!(v, v.tune.logf; args...)

function sample!(v::AMMVariate, logf::Function; adapt::Bool=true)
  tune = v.tune
  n = length(v)

  setadapt!(v, adapt)

  x = tune.SigmaL * randn(n)
  if tune.m > 2 * n
    x = tune.beta * x + (1.0 - tune.beta) * (tune.SigmaLm * randn(n))
  end
  x += v
  if rand() < exp(logf(x) - logf(v.value))
    v[:] = x
  end

  if tune.adapt
    tune.m += 1
    p = tune.m / (tune.m + 1.0)
    tune.Mv = p * tune.Mv + (1.0 - p) * v
    tune.Mvv = p * tune.Mvv + (1.0 - p) * v * v'
    Sigma = (tune.scale^2 / n / p) * (tune.Mvv - tune.Mv * tune.Mv')
    F = cholfact(Hermitian(Sigma), Val{true})
    if rank(F) == n
      tune.SigmaLm = F[:P] * F[:L]
    end
  end

  v
end


function setadapt!(v::AMMVariate, adapt::Bool)
  tune = v.tune
  if adapt && !tune.adapt
    n = length(v)
    tune.m = 0
    tune.Mv = v
    tune.Mvv = v * v'
    tune.SigmaLm = zeros(n, n)
  end
  tune.adapt = adapt
  v
end
