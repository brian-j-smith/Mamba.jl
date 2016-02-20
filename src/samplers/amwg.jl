#################### Adaptive Metropolis within Gibbs ####################

#################### Types and Constructors ####################

type AMWGTune <: SamplerTune
  logf::Nullable{Function}
  adapt::Bool
  accept::Vector{Int}
  batchsize::Int
  m::Int
  sigma::Vector{Float64}
  target::Float64

  AMWGTune() = new()

  function AMWGTune{T<:Real}(x::Vector, sigma::Vector{T},
                             logf::Nullable{Function}; batchsize::Integer=50,
                             target::Real=0.44)
    new(logf, false, zeros(Int, length(x)), batchsize, 0, copy(sigma), target)
  end
end

AMWGTune{T<:Real}(x::Vector, sigma::ElementOrVector{T}; args...) =
  AMWGTune(x, sigma, Nullable{Function}(); args...)

function AMWGTune{T<:Real}(x::Vector, sigma::ElementOrVector{T}, logf::Function;
                           args...)
  AMWGTune(x, sigma, Nullable{Function}(logf); args...)
end

AMWGTune(x::Vector, sigma::Real, logf::Nullable{Function}; args...) =
  AMWGTune(x, fill(sigma, length(x)), logf; args...)


typealias AMWGVariate SamplerVariate{AMWGTune}

function validate(v::AMWGVariate)
  n = length(v)
  length(v.tune.sigma) == n ||
    throw(ArgumentError("length(sigma) differs from variate length $n"))
  v
end


#################### Sampler Constructor ####################

function AMWG{T<:Real}(params::ElementOrVector{Symbol},
                       sigma::ElementOrVector{T}; adapt::Symbol=:all, args...)
  adapt in [:all, :burnin, :none] ||
    throw(ArgumentError("adapt must be one of :all, :burnin, or :none"))

  samplerfx = function(model::Model, block::Integer)
    block = SamplingBlock(model, block, true)
    v = SamplerVariate(block, sigma; args...)
    isadapt = adapt == :burnin ? model.iter <= model.burnin :
              adapt == :all ? true : false
    sample!(v, x -> logpdf!(block, x), adapt=isadapt)
    relist(block, v)
  end
  Sampler(params, samplerfx, AMWGTune())
end


#################### Sampling Functions ####################

sample!(v::AMWGVariate; args...) = sample!(v, v.tune.logf; args...)

function sample!(v::AMWGVariate, logf::Function; adapt::Bool=true)
  tune = v.tune
  setadapt!(v, adapt)
  if tune.adapt
    tune.m += 1
    amwg_sub!(v, logf)
    if tune.m % tune.batchsize == 0
      delta = min(0.01, (tune.m / tune.batchsize)^-0.5)
      for i in 1:length(tune.sigma)
        epsilon = tune.accept[i] / tune.m < tune.target ? -delta : delta
        tune.sigma[i] *= exp(epsilon)
      end
    end
  else
    amwg_sub!(v, logf)
  end
  v
end


function setadapt!(v::AMWGVariate, adapt::Bool)
  tune = v.tune
  if adapt && !tune.adapt
    tune.accept[:] = 0
    tune.m = 0
  end
  tune.adapt = adapt
  v
end


function amwg_sub!(v::AMWGVariate, logf::Function)
  logf0 = logf(v.value)
  n = length(v)
  z = v.tune.sigma .* randn(n)
  for i in 1:n
    x = v[i]
    v[i] += z[i]
    logfprime = logf(v.value)
    if rand() < exp(logfprime - logf0)
      logf0 = logfprime
      v.tune.accept[i] += v.tune.adapt
    else
      v[i] = x
    end
  end
  v
end
