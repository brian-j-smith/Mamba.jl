#################### No-U-Turn Sampler ####################

#################### Types and Constructors ####################

type NUTSTune <: SamplerTune
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

  function NUTSTune(value::Vector{Float64}=Float64[])
    new(
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
  end
end


typealias NUTSVariate SamplerVariate{NUTSTune}


#################### Sampler Constructor ####################

function NUTS(params::ElementOrVector{Symbol}; dtype::Symbol=:forward,
              target::Real=0.6)
  samplerfx = function(model::Model, block::Integer)
    v = SamplerVariate(model, block, true)
    tune = v.tune
    if model.iter == 1
      f = x -> logpdfgrad(model, x, block, dtype)
      tune.epsilon = nutsepsilon(v, f)
    end
    f = x -> logpdfgrad!(model, x, block, dtype)
    nuts!(v, tune.epsilon, f, adapt=model.iter <= model.burnin, target=target)
    relist(model, v, block, true)
  end
  Sampler(params, samplerfx, NUTSTune())
end


#################### Sampling Functions ####################

function nutsepsilon(v::NUTSVariate, logfgrad::Function)
  n = length(v)
  _, r0, logf0, grad0 = leapfrog(v.value, randn(n), zeros(n), 0.0, logfgrad)
  epsilon = 1.0
  _, rprime, logfprime, gradprime = leapfrog(v.value, r0, grad0, epsilon,
                                             logfgrad)
  prob = exp(logfprime - logf0 - 0.5 * (dot(rprime) - dot(r0)))
  pm = 2 * (prob > 0.5) - 1
  while prob^pm > 0.5^pm
    epsilon *= 2.0^pm
    _, rprime, logfprime, _ = leapfrog(v.value, r0, grad0, epsilon, logfgrad)
    prob = exp(logfprime - logf0 - 0.5 * (dot(rprime) - dot(r0)))
  end
  epsilon
end


function nuts!(v::NUTSVariate, epsilon::Real, logfgrad::Function;
               adapt::Bool=false, target::Real=0.6)
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
    nuts_sub!(v, tune.epsilon, logfgrad)
    p = 1.0 / (tune.m + tune.t0)
    tune.Hbar = (1.0 - p) * tune.Hbar +
                p * (tune.target - tune.alpha / tune.nalpha)
    tune.epsilon = exp(tune.mu - sqrt(tune.m) * tune.Hbar / tune.gamma)
    p = tune.m^-tune.kappa
    tune.epsilonbar = exp(p * log(tune.epsilon) +
                          (1.0 - p) * log(tune.epsilonbar))
  else
    tune.epsilon = tune.adapt ? tune.epsilonbar : epsilon
    nuts_sub!(v, tune.epsilon, logfgrad)
  end

  v
end


function nuts_sub!(v::NUTSVariate, epsilon::Real, logfgrad::Function)
  n = length(v)
  x, r, logf, grad = leapfrog(v.value, randn(n), zeros(n), 0.0, logfgrad)
  logp0 = logf - 0.5 * dot(r)
  logu0 = logp0 + log(rand())
  xminus = xplus = x
  rminus = rplus = r
  gradminus = gradplus = grad
  j = 0
  n = 1
  s = true
  while s
    pm = 2 * (rand() > 0.5) - 1
    if pm == -1
      xminus, rminus, gradminus, _, _, _, xprime, nprime, sprime, alpha,
        nalpha = buildtree(xminus, rminus, gradminus, pm, j, epsilon, logfgrad,
                           logp0, logu0)
    else
      _, _, _, xplus, rplus, gradplus, xprime, nprime, sprime, alpha, nalpha =
        buildtree(xplus, rplus, gradplus, pm, j, epsilon, logfgrad, logp0,
                  logu0)
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
                  epsilon::Real, logfgrad::Function)
  r += (0.5 * epsilon) * grad
  x += epsilon * r
  logf, grad = logfgrad(x)
  r += (0.5 * epsilon) * grad
  x, r, logf, grad
end


function buildtree(x::Vector{Float64}, r::Vector{Float64},
                   grad::Vector{Float64}, pm::Integer, j::Integer,
                   epsilon::Real, logfgrad::Function, logp0::Real, logu0::Real)
  if j == 0
    xprime, rprime, logfprime, gradprime = leapfrog(x, r, grad, pm * epsilon,
                                                    logfgrad)
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
      alphaprime, nalphaprime = buildtree(x, r, grad, pm, j - 1, epsilon,
                                          logfgrad, logp0, logu0)
    if sprime
      if pm == -1
        xminus, rminus, gradminus, _, _, _, xprime2, nprime2, sprime2,
          alphaprime2, nalphaprime2 = buildtree(xminus, rminus, gradminus, pm,
                                                j - 1, epsilon, logfgrad, logp0,
                                                logu0)
      else
        _, _, _, xplus, rplus, gradplus, xprime2, nprime2, sprime2,
          alphaprime2, nalphaprime2 = buildtree(xplus, rplus, gradplus, pm,
                                                j - 1, epsilon, logfgrad, logp0,
                                                logu0)
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


function nouturn(xminus::Vector{Float64}, xplus::Vector{Float64},
                 rminus::Vector{Float64}, rplus::Vector{Float64})
  xdiff = xplus - xminus
  dot(xdiff, rminus) >= 0 && dot(xdiff, rplus) >= 0
end
