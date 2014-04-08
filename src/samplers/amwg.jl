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

  function VariateAMWG{T<:Real}(x::Vector{T}, tune::TuneAMWG)
    new(VariateType[x...], tune)
  end
end

function VariateAMWG{T<:Real}(x::Vector{T}, tune=nothing)
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


#################### Sampling Functions ####################

function amwg{T<:Real}(x::Vector{T}, sigma::Vector{Float64}, logf::Function;
           adapt::Bool=false, batch::Integer=50, target::Real=0.44)
  amwg!(VariateAMWG(x), sigma, logf, adapt=adapt, batch=batch,
        target=target)
end

function amwg!(v::VariateAMWG, sigma::Vector{Float64}, logf::Function;
           adapt::Bool=false, batch::Integer=50, target::Real=0.44)
  tune = v.tune

  if adapt
    if !tune.adapt
      tune.adapt = true
      tune.accept[:] = 0
      tune.batch = batch
      tune.m = 0
      tune.sigma = sigma
      tune.target = target
    end
    tune.m += 1
    amwg_sub!(v, tune.sigma, logf)
    if tune.m % tune.batch == 0
      delta = min(0.01, (tune.m / tune.batch)^-0.5)
      for i in 1:length(tune.sigma)
        tune.sigma[i] *= tune.accept[i] / tune.batch < tune.target ?
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
  logf0 = logf(v.data)
  d = length(v)
  z = randn(d) .* sigma
  for i in 1:d
    x = v[i]
    v[i] += z[i]
    logfprime = logf(v.data)
    if rand() < exp(logfprime - logf0)
      logf0 = logfprime
      v.tune.accept[i] += 1
    else
      v[i] = x
    end
  end
  v
end
