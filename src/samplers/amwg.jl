#################### Adaptive Metropolis within Gibbs ####################

#################### Types and Constructors ####################

type AMWGTune <: SamplerTune
  adapt::Bool
  accept::Vector{Int}
  batchsize::Int
  m::Int
  sigma::Vector{Float64}
  target::Real
end

type AMWGVariate <: SamplerVariate
  value::Vector{Float64}
  tune::AMWGTune

  AMWGVariate{T<:Real}(x::AbstractVector{T}, tune::AMWGTune) = new(x, tune)
end

function AMWGVariate{T<:Real}(x::AbstractVector{T}, tune=nothing)
  tune = AMWGTune(
    false,
    zeros(Int, length(x)),
    50,
    0,
    Array(Float64, 0),
    0.44
  )
  AMWGVariate(x, tune)
end


#################### Sampler Constructor ####################

function AMWG{T<:Real}(params::Vector{Symbol}, sigma::Vector{T};
                       adapt::Symbol=:all, batchsize::Integer=50,
                       target::Real=0.44)
  adapt in [:all, :burnin, :none] ||
    throw(ArgumentError("adapt must be one of :all, :burnin, or :none"))

  sigma = Float64[sigma...]
  samplerfx = function(model::Model, block::Integer)
    tunepar = tune(model, block)
    x = unlist(model, block, true)
    v = AMWGVariate(x, tunepar["sampler"])
    f = x -> logpdf!(model, x, block, true)
    isadapt = adapt == :burnin ? model.iter <= model.burnin :
              adapt == :all ? true : false
    amwg!(v, sigma, f, adapt=isadapt, batchsize=batchsize, target=target)
    tunepar["sampler"] = v.tune
    relist(model, v, block, true)
  end
  Sampler(params, samplerfx, Dict{AbstractString, Any}("sampler" => nothing))
end


#################### Sampling Functions ####################

function amwg!(v::AMWGVariate, sigma::Vector{Float64}, logf::Function;
               adapt::Bool=true, batchsize::Integer=50, target::Real=0.44)
  tune = v.tune

  if adapt
    if !tune.adapt
      tune.adapt = true
      tune.accept[:] = 0
      tune.batchsize = batchsize
      tune.m = 0
      tune.sigma = copy(sigma)
      tune.target = target
    end
    tune.m += 1
    amwg_sub!(v, tune.sigma, logf)
    if tune.m % tune.batchsize == 0
      delta = min(0.01, (tune.m / tune.batchsize)^-0.5)
      for i in 1:length(tune.sigma)
        epsilon = tune.accept[i] / tune.m < tune.target ? -delta : delta
        tune.sigma[i] *= exp(epsilon)
      end
    end
  else
    if !tune.adapt
      tune.sigma = sigma
    end
    amwg_sub!(v, tune.sigma, logf)
  end

  v
end

function amwg_sub!(v::AMWGVariate, sigma::Vector{Float64}, logf::Function)
  logf0 = logf(v.value)
  d = length(v)
  z = randn(d) .* sigma
  for i in 1:d
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
