#################### Adaptive Metropolis within Gibbs ####################

#################### Types ####################

type AMWGTune
  adapt::Bool
  accept::Vector{Integer}
  batchsize::Integer
  m::Integer
  sigma::Vector{Float64}
  target::Real
end

type AMWGVariate <: VectorVariate
  value::Vector{VariateType}
  tune::AMWGTune

  AMWGVariate(x::Vector{VariateType}, tune::AMWGTune) = new(x, tune)
end

function AMWGVariate(x::Vector{VariateType}, tune=nothing)
  tune = AMWGTune(
    false,
    zeros(Integer, length(x)),
    50,
    0,
    Array(Float64, 0),
    0.44
  )
  AMWGVariate(x, tune)
end


#################### Sampler Constructor ####################

function AMWG{T<:Real}(params::Vector{Symbol}, sigma::Vector{T};
           adapt::Symbol=:all, batchsize::Integer=50, target::Real=0.44)
  in(adapt, [:all, :burnin, :none]) ||
    error("adapt argument must be one of :all, :burnin, or :none")

  Sampler(params,
    quote
      x = unlist(model, block, true)
      tunepar = tune(model, block)
      v = AMWGVariate(x, tunepar["sampler"])
      adapt = tunepar["adapt"] == :burnin ? model.iter <= model.burnin :
              tunepar["adapt"] == :all ? true : false
      f = x -> logpdf!(model, x, block, true)
      amwg!(v, tunepar["sigma"], f, adapt=adapt, batchsize=tunepar["batchsize"],
            target=tunepar["target"])
      tunepar["sampler"] = v.tune
      relist(model, v.value, block, true)
    end,
    ["sigma" => Float64[sigma...], "adapt" => adapt, "batchsize" => batchsize,
     "target" => target, "sampler" => nothing]
  )
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
      tune.sigma = sigma
      tune.target = target
    end
    tune.m += 1
    amwg_sub!(v, tune.sigma, logf)
    if tune.m % tune.batchsize == 0
      delta = min(0.01, (tune.m / tune.batchsize)^-0.5)
      for i in 1:length(tune.sigma)
        tune.sigma[i] *= tune.accept[i] / tune.batchsize < tune.target ?
          exp(-delta) : exp(delta)
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
      v.tune.accept[i] += 1
    else
      v[i] = x
    end
  end
  v
end
