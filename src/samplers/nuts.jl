#################### No-U-Turn Sampler ####################

#################### Types and Constructors ####################

type NUTSTune
  adapt::Bool
  alpha::Float64
  epsilon::Float64
  epsilonbar::Float64
  gamma::Float64
  Hbar::Float64
  kappa::Float64
  m::Int
  mu::Float64
  nalpha::Int
  t0::Float64
  target::Float64
end

type NUTSVariate <: VectorVariate
  value::Vector{Float64}
  tune::NUTSTune

  NUTSVariate(x::Vector{Float64}, tune::NUTSTune) = new(x, tune)
end

function NUTSVariate(x::Vector{Float64}, tune=nothing)
  tune = NUTSTune(
    false,
    0.0,
    NaN,
    1.0,
    0.05,
    0.0,
    0.75,
    0,
    NaN,
    0,
    10.0,
    0.44
  )
  NUTSVariate(x, tune)
end


#################### Sampler Constructor ####################

function NUTS(params::Vector{Symbol}; dtype::Symbol=:forward, target::Real=0.6)
  Sampler(params,
    quote
      tunepar = tune(model, block)
      x = unlist(model, block, true)
      v = NUTSVariate(x, tunepar["sampler"])
      if model.iter <= 1
        f = x -> nutsfx(model, x, block, tunepar["dtype"])
        tunepar["epsilon"] = nutsepsilon(v, f)
      end
      f = x -> nutsfx!(model, x, block, tunepar["dtype"])
      nuts!(v, tunepar["epsilon"], f, adapt=model.iter <= model.burnin,
            target=tunepar["target"])
      tunepar["sampler"] = v.tune
      relist(model, v, block, true)
    end,
    Dict("epsilon" => 1.0, "target" => target, "dtype" => dtype,
         "sampler" => nothing)
  )
end

function nutsfx{T<:Real}(m::Model, x::Vector{T}, block::Integer, dtype::Symbol)
  logf = logpdf(m, x, block, true)
  grad = isfinite(logf) ?
           gradlogpdf(m, x, block, true, dtype=dtype) :
           zeros(length(x))
  logf, grad
end

function nutsfx!{T<:Real}(m::Model, x::Vector{T}, block::Integer, dtype::Symbol)
  logf = logpdf!(m, x, block, true)
  grad = isfinite(logf) ?
           gradlogpdf!(m, x, block, true, dtype=dtype) :
           zeros(length(x))
  logf, grad
end

export nutsfx, nutsfx!


#################### Sampling Functions ####################

function nutsepsilon(v::NUTSVariate, fx::Function)
  n = length(v)
  _, r0, grad0, logf0 = leapfrog(v.value, randn(n), zeros(n), 0.0, fx)
  epsilon = 1.0
  _, rprime, gradprime, logfprime = leapfrog(v.value, r0, grad0, epsilon, fx)
  prob = exp(logfprime - logf0 - 0.5 * (dot(rprime) - dot(r0)))
  d = 2 * (prob > 0.5) - 1
  while prob^d > 0.5^d
    epsilon *= 2.0^d
    _, rprime, _, logfprime = leapfrog(v.value, r0, grad0, epsilon, fx)
    prob = exp(logfprime - logf0 - 0.5 * (dot(rprime) - dot(r0)))
  end
  epsilon
end

function nuts!(v::NUTSVariate, epsilon::Real, fx::Function; adapt::Bool=false,
               target::Real=0.6)
  tune = v.tune

  if adapt
    if !tune.adapt
      tune.adapt = true
      tune.epsilon = epsilon
      tune.m = 0
      tune.mu = log(10.0 * epsilon)
      tune.target = target
    end
    tune.m += 1
    nuts_sub!(v, tune.epsilon, fx)
    p = 1.0 / (tune.m + tune.t0)
    tune.Hbar = (1.0 - p) * tune.Hbar +
                p * (tune.target - tune.alpha / tune.nalpha)
    tune.epsilon = exp(tune.mu - sqrt(tune.m) * tune.Hbar / tune.gamma)
    p = tune.m^-tune.kappa
    tune.epsilonbar = exp(p * log(tune.epsilon) +
                          (1.0 - p) * log(tune.epsilonbar))
  else
    tune.epsilon = tune.adapt ? tune.epsilonbar : epsilon
    nuts_sub!(v, tune.epsilon, fx)
  end

  v
end

function nuts_sub!(v::NUTSVariate, epsilon::Real, fx::Function)
  n = length(v)
  x, r, grad, logf = leapfrog(v.value, randn(n), zeros(n), 0.0, fx)
  logp0 = logf - 0.5 * dot(r)
  logu0 = logp0 + log(rand())
  xminus = xplus = x
  rminus = rplus = r
  gradminus = gradplus = grad
  j = 0
  n = 1
  s = true
  while s
    d = 2 * (rand() > 0.5) - 1
    if d == -1
      xminus, rminus, gradminus, _, _, _, xprime, nprime, sprime, alpha, nalpha =
        buildtree(xminus, rminus, gradminus, d, j, epsilon, fx, logp0, logu0)
    else
      _, _, _, xplus, rplus, gradplus, xprime, nprime, sprime, alpha, nalpha =
        buildtree(xplus, rplus, gradplus, d, j, epsilon, fx, logp0, logu0)
    end
    if sprime && rand() < nprime / n
      v[:] = xprime
    end
    j += 1
    n += nprime
    s = sprime && nouturn(xminus, xplus, rminus, rplus)
    v.tune.alpha, v.tune.nalpha = alpha, nalpha
  end
  v
end

function leapfrog(x::Vector{Float64}, r::Vector{Float64}, grad::Vector{Float64},
                  epsilon::Real, fx::Function)
  r += (0.5 * epsilon) * grad
  x += epsilon * r
  logf, grad = fx(x)
  r += (0.5 * epsilon) * grad
  x, r, grad, logf
end

function buildtree(x::Vector, r::Vector, grad::Vector, d::Integer, j::Integer,
                   epsilon::Real, fx::Function, logp0::Real, logu0::Real)
  if j == 0
    xprime, rprime, gradprime, logfprime = leapfrog(x, r, grad, d * epsilon, fx)
    logpprime = logfprime - 0.5 * dot(rprime)
    nprime = Int(logu0 < logpprime)
    sprime = logu0 < logpprime + 1000.0
    xminus = xplus = xprime
    rminus = rplus = rprime
    gradminus = gradplus = gradprime
    alphaprime = min(1.0, exp(logpprime - logp0))
    nalphaprime = 1
  else
    xminus, rminus, gradminus, xplus, rplus, gradplus, xprime, nprime, sprime,
      alphaprime, nalphaprime =
        buildtree(x, r, grad, d, j-1, epsilon, fx, logp0, logu0)
    if sprime
      if d == -1
        xminus, rminus, gradminus, _, _, _, xprime2, nprime2, sprime2,
          alphaprime2, nalphaprime2 =
            buildtree(xminus, rminus, gradminus, d, j-1, epsilon, fx, logp0, logu0)
      else
        _, _, _, xplus, rplus, gradplus, xprime2, nprime2, sprime2,
          alphaprime2, nalphaprime2 =
            buildtree(xplus, rplus, gradplus, d, j-1, epsilon, fx, logp0, logu0)
      end
      if rand() < nprime2 / (nprime + nprime2)
        xprime = xprime2
      end
      nprime += nprime2
      sprime = sprime2 && nouturn(xminus, xplus, rminus, rplus)
      alphaprime += alphaprime2
      nalphaprime += nalphaprime2
    end
  end
  xminus, rminus, gradminus, xplus, rplus, gradplus, xprime, nprime, sprime,
    alphaprime, nalphaprime
end

function nouturn(xminus, xplus, rminus, rplus)
  xdiff = xplus - xminus
  dot(xdiff, rminus) >= 0 && dot(xdiff, rplus) >= 0
end
