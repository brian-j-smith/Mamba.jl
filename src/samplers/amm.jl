#################### Adaptive Mixture Metropolis ####################

#################### Types ####################

type TuneAMM
  adapt::Bool
  beta::Real
  m::Integer
  mu::Vector{Float64}
  scale::Real
  Sigma::Cholesky{Float64}
  SigmaLm::Matrix{Float64}
end

type VariateAMM <: VariateVector
  data::Vector{VariateType}
  tune::TuneAMM

  function VariateAMM(x::Vector, tune::TuneAMM)
    new(VariateType[x...], tune)
  end
end

function VariateAMM(x::Vector, Sigma::Cholesky)
  n = length(x)
  size(Sigma)[1] == n || error("x and Sigma dimensions must match")
  tune = TuneAMM(
    false,
    0.05,
    0,
    zeros(n),
    2.38^2,
    Sigma,
    zeros(n,n)
  )
  VariateAMM(x, tune)
end

function VariateAMM(x::Vector, tune=nothing)
  VariateAMM(x, Cholesky(eye(length(x)), 'U'))
end


#################### Sampling Functions ####################

function amm(x::Vector, Sigma::Union(Matrix, Cholesky), logf::Function, args...;
             adapt::Bool=false)
  amm!(VariateAMM(x), Sigma, logf, args...)
end

function amm!(v::VariateAMM, Sigma::Matrix, logf::Function, args...;
              adapt::Bool=false)
  amm!(v, cholfact(Sigma), logf, args...)
end

function amm!(v::VariateAMM, Sigma::Cholesky, logf::Function, args...;
              adapt::Bool=false)
  tune = v.tune

  d = length(v)
  if adapt
    if !tune.adapt
      tune.adapt = true
      tune.m = 0
      tune.mu[:] = 0
      tune.Sigma = Sigma
      tune.SigmaLm = zeros(d,d)
    end
    tune.m += 1
    x = v + tune.Sigma[:L] * randn(d)
    if tune.m > 2 * d
      x = tune.beta * x + (1 - tune.beta) * (v + tune.SigmaLm * randn(d))
    end
    if rand() < exp(logf(x, args...) - logf(v.data, args...))
      v[:] = x
    end
    sd = tune.scale / d
    p = 1 / (tune.m + 1)
    mu = (1 - p) * tune.mu + p * v
    Sigma = (1 - 1 / tune.m) * tune.SigmaLm * tune.SigmaLm' +
            sd * tune.mu * tune.mu' - sd * (1 + 1 / tune.m) * mu * mu' +
            sd / tune.m * v * v'
    tune.mu[:] = mu
    F = cholpfact(Sigma)
    tune.SigmaLm = F[:P] * F[:L]
  else
    if tune.adapt
      x = v + tune.beta * (tune.Sigma[:L] * randn(d)) +
          (1 - tune.beta) * (tune.SigmaLm * randn(d))
    else
      tune.Sigma = Sigma
      x = v + tune.Sigma[:L] * randn(d)
    end
    if rand() < exp(logf(x, args...) - logf(v.data, args...))
      v[:] = x
    end
  end

  v
end
