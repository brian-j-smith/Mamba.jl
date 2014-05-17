#################### Adaptive Metropolis within Gibbs ####################

#################### Types ####################

type TuneAMWG
  adapt::Bool
  accept::Vector{Integer}
  batchsize::Integer
  m::Integer
  sigma::Vector{Float64}
  target::Real
end

type VariateAMWG <: VariateVector
  value::Vector{VariateType}
  tune::TuneAMWG

  VariateAMWG(x::Vector{VariateType}, tune::TuneAMWG) = new(x, tune)
end

function VariateAMWG(x::Vector{VariateType}, tune=nothing)
  tune = TuneAMWG(
    false,
    zeros(Integer, length(x)),
    50,
    0,
    Array(Float64, 0),
    0.44
  )
  VariateAMWG(x, tune)
end


#################### MCMCSampler Constructor ####################

function AMWG{T<:String,U<:Real}(params::Vector{T}, sigma::Vector{U};
           adapt::Symbol=:all, batchsize::Integer=50, target::Real=0.44)
  any(adapt .== [:all, :burnin, :none]) ||
    error("adapt argument must be one of :all, :burnin, or :none")

  MCMCSampler(params,
    quote
      x = unlist(model, block, true)
      tunepar = tune(model, block)
      v = VariateAMWG(x, tunepar["sampler"])
      adapt = tunepar["adapt"] == :burnin ? model.iter <= model.burnin :
              tunepar["adapt"] == :all ? true : false
      f = x -> logpdf!(model, x, block, true)
      amwg!(v, tunepar["sigma"], f, adapt=adapt, batchsize=tunepar["batchsize"],
            target=tunepar["target"])
      tunepar["sampler"] = v.tune
      relist(model, v.value, block, true)
    end,
    ["sigma" => sigma, "adapt" => adapt, "batchsize" => batchsize,
     "target" => target, "sampler" => nothing]
  )
end


#################### Sampling Functions ####################

function amwg!(v::VariateAMWG, sigma::Vector{Float64}, logf::Function;
           adapt::Bool=false, batchsize::Integer=50, target::Real=0.44)
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

function amwg_sub!(v::VariateAMWG, sigma::Vector{Float64}, logf::Function)
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
