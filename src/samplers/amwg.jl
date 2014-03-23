#################### Adaptive Metropolis within Gibbs ####################

#################### Types ####################

type TuneAMWG
  adapt::Bool
  accept::Vector{Integer}
  batch::Integer
  m::Integer
  sigma::Vector{Float64}
  target::Real
end

type VariateAMWG <: VariateVector
  data::Vector{VariateType}
  tune::TuneAMWG

  function VariateAMWG(x::Vector, tune::TuneAMWG)
    new(VariateType[x...], tune)
  end
end

function VariateAMWG(x::Vector, sigma::Vector; batch::Integer=50,
                     target::Real=0.44)
  n = length(x)
  length(sigma) == n || error("x and sigma dimensions must match")
  tune = TuneAMWG(
    false,
    Integer[zeros(n)...],
    batch,
    0,
    Float64[sigma...],
    target
  )
  VariateAMWG(x, tune)
end

function VariateAMWG(x::Vector, tune=nothing)
  VariateAMWG(x, ones(length(x)))
end


#################### Sampling Functions ####################

function amwg(x::Vector, sigma::Vector, logf::Function, args...;
              adapt::Bool=false, batch::Integer=50, target::Real=0.44)
  amwg!(VariateAMWG(x), sigma, logf, args..., adapt=adapt, batch=batch,
        target=target)
end

function amwg!(v::VariateAMWG, sigma::Vector, logf::Function, args...;
               adapt::Bool=false, batch::Integer=50, target::Real=0.44)
  tune = v.tune

  if adapt
    if !tune.adapt
      tune.adapt = true
      tune.accept[:] = 0
      tune.batch = batch
      tune.m = 0
      tune.sigma[:] = sigma
      tune.target = target
    end
    tune.m += 1
    amwg_sub!(v, tune.sigma, logf, args...)
    if tune.m % tune.batch == 0
      delta = min(0.01, (tune.m / tune.batch)^-0.5)
      for i in 1:length(tune.sigma)
        tune.sigma[i] *= tune.accept[i] / tune.batch < tune.target ?
          exp(-delta) : exp(delta)
      end
    end
  else
    if !tune.adapt
      tune.sigma[:] = sigma
    end
    amwg_sub!(v, tune.sigma, logf, args...)
  end

  v
end

function amwg_sub!(v::VariateAMWG, sigma::Vector, logf::Function, args...)
  logf0 = logf(v.data, args...)
  d = length(v)
  z = randn(d) .* sigma
  for i in 1:d
    x = v[i]
    v[i] += z[i]
    logfprime = logf(v.data, args...)
    if rand() < exp(logfprime - logf0)
      logf0 = logfprime
      v.tune.accept[i] += 1
    else
      v[i] = x
    end
  end
  v
end
