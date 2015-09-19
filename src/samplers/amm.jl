#################### Adaptive Mixture Metropolis ####################

#################### Types ####################

type AMMTune
  adapt::Bool
  beta::Real
  m::Int
  Mv::Vector{Float64}
  Mvv::Matrix{Float64}
  scale::Real
  SigmaF::Cholesky{Float64}
  SigmaLm::Matrix{Float64}
end

type AMMVariate <: VectorVariate
  value::Vector{Float64}
  tune::AMMTune

  AMMVariate(x::Vector{Float64}, tune::AMMTune) = new(x, tune)
end

function AMMVariate(x::Vector{Float64}, tune=nothing)
  tune = AMMTune(
    false,
    0.05,
    0,
    Array(Float64, 0),
    Array(Float64, 0, 0),
    2.38^2,
    Cholesky(Array(Float64, 0, 0), :U),
    Array(Float64, 0, 0)
  )
  AMMVariate(x, tune)
end


#################### Sampler Constructor ####################

function AMM{T<:Real}(params::Vector{Symbol}, Sigma::Matrix{T};
                      adapt::Symbol=:all)
  in(adapt, [:all, :burnin, :none]) ||
    error("adapt argument must be one of :all, :burnin, or :none")

  Sampler(params,
    quote
      tunepar = tune(model, block)
      x = unlist(model, block, true)
      v = AMMVariate(x, tunepar["sampler"])
      adapt = tunepar["adapt"] == :burnin ? model.iter <= model.burnin :
              tunepar["adapt"] == :all ? true : false
      f = x -> logpdf!(model, x, block, true)
      amm!(v, tunepar["SigmaF"], f, adapt=adapt)
      tunepar["sampler"] = v.tune
      relist(model, v.value, block, true)
    end,
    Dict("SigmaF" => cholfact(Sigma), "adapt" => adapt, "sampler" => nothing)
  )
end


#################### Sampling Functions ####################

function amm!(v::AMMVariate, SigmaF::Cholesky{Float64}, logf::Function;
              adapt::Bool=true)
  tune = v.tune

  d = length(v)
  if adapt
    if !tune.adapt
      tune.adapt = true
      tune.m = 0
      tune.Mv = v.value
      tune.Mvv = v * v'
      tune.SigmaF = SigmaF
      tune.SigmaLm = zeros(d, d)
    end
    x = v + tune.SigmaF[:L] * randn(d)
    if tune.m > 2 * d
      x = tune.beta * x + (1.0 - tune.beta) * (v + tune.SigmaLm * randn(d))
    end
    if rand() < exp(logf(x) - logf(v.value))
      v[:] = x
    end
    tune.m += 1
    sd = tune.scale / d
    p = tune.m / (tune.m + 1.0)
    tune.Mv = p * tune.Mv + (1.0 - p) * v
    tune.Mvv = p * tune.Mvv + (1.0 - p) * v * v'
    Sigma = (sd / p) * (tune.Mvv - tune.Mv * tune.Mv')
    F = cholfact(Sigma, :U, Val{true})
    if rank(F) == d
      tune.SigmaLm = F[:P] * F[:L]
    end
  else
    if tune.adapt
      x = v + tune.beta * (tune.SigmaF[:L] * randn(d)) +
          (1.0 - tune.beta) * (tune.SigmaLm * randn(d))
    else
      tune.SigmaF = SigmaF
      x = v + tune.SigmaF[:L] * randn(d)
    end
    if rand() < exp(logf(x) - logf(v.value))
      v[:] = x
    end
  end

  v
end
