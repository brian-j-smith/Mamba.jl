#################### No-U-Turn Sampler ####################

#################### Types ####################

type TuneNUTS
  adapt::Bool
  alpha::Float64
  eps::Float64
  epsbar::Float64
  gamma::Float64
  Hbar::Float64
  kappa::Float64
  m::Integer
  mu::Float64
  nalpha::Integer
  t0::Float64
  target::Float64
end

type VariateNUTS <: VariateVector
  value::Vector{VariateType}
  tune::TuneNUTS
end

function VariateNUTS(x::Vector{VariateType}, tune=nothing)
  tune = TuneNUTS(
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
  VariateNUTS(x, tune)
end


#################### MCMCSampler Constructor ####################

function NUTS{T<:String}(params::Vector{T}; dtype::Symbol=:forward,
           target::Real=0.6)
  MCMCSampler(params,
    quote
      x = unlist(model, block, true)
      tunepar = tune(model, block)
      v = VariateNUTS(x, tunepar["sampler"])
      if model.iter <= 1
        f = x -> nutsfx(model, x, block, true, tunepar["dtype"])
        tunepar["eps"] = nutseps(v, f)
      end
      f = x -> nutsfx!(model, x, block, true, tunepar["dtype"])
      nuts!(v, tunepar["eps"], f, adapt=model.iter <= model.burnin,
            target=tunepar["target"])
      tunepar["sampler"] = v.tune
      relist(model, v.value, block, true)
    end,
    ["eps" => 1.0, "target" => target, "dtype" => dtype, "sampler" => nothing]
  )
end

function nutsfx{T<:Real}(m::MCMCModel, x::Vector{T}, block::Integer,
           transform::Bool, dtype::Symbol)
  a = logpdf(m, x, block, transform)
  b = gradient(m, x, block, transform, dtype)
  a, b
end

function nutsfx!{T<:Real}(m::MCMCModel, x::Vector{T}, block::Integer,
           transform::Bool, dtype::Symbol)
  a = logpdf!(m, x, block, transform)
  b = gradient(m, block, transform, dtype)
  a, b
end

export nutsfx, nutsfx!


#################### Sampling Functions ####################

function nutseps(v::VariateNUTS, fx::Function)
  d = length(v)
  node0 = leapfrog(v.value, randn(d), 0.0, zeros(d), fx)
  eps = 1.0
  node = leapfrog(v.value, node0[:r], eps, node0[:grad], fx)
  p = exp(node[:logf] - node0[:logf] - 0.5 * (dot(node[:r]) - dot(node0[:r])))
  a = p > 0.5 ? 1 : -1
  while p^a > 2.0^-a
    eps *= 2.0^a
    node = leapfrog(v.value, node0[:r], eps, node0[:grad], fx)
    p = exp(node[:logf] - node0[:logf] - 0.5 * (dot(node[:r]) - dot(node0[:r])))
  end
  eps
end

function leapfrog{T<:Real,U<:Real,V<:Real}(x::Vector{T}, r::Vector{U},
           eps::Real, grad::Vector{V}, fx::Function)
  r += (eps / 2.0) * grad
  x += eps * r
  logf, grad = fx(x)
  r += (eps / 2.0) * grad
  [:x => x, :r => r, :logf => logf, :grad => grad]
end

function nuts!(v::VariateNUTS, eps::Real, fx::Function; adapt::Bool=false,
           target::Real=0.6)
  tune = v.tune

  if adapt
    if !tune.adapt
      tune.adapt = true
      tune.eps = eps
      tune.m = 0
      tune.mu = log(10.0 * eps)
      tune.target = target
    end
    tune.m += 1
    nuts_sub!(v, tune.eps, fx)
    p = 1.0 / (tune.m + tune.t0)
    tune.Hbar = (1 - p) * tune.Hbar +
                p * (tune.target - tune.alpha / tune.nalpha)
    tune.eps = exp(tune.mu - sqrt(tune.m) * tune.Hbar / tune.gamma)
    p = tune.m^-tune.kappa
    tune.epsbar = exp(p * log(tune.eps) + (1 - p) * log(tune.epsbar))
  else
    tune.eps = tune.adapt ? tune.epsbar : eps
    nuts_sub!(v, tune.eps, fx)
  end

  v
end

function nuts_sub!(v::VariateNUTS, eps::Real, fx::Function)
  d = length(v)
  node0 = leapfrog(v.value, randn(d), 0.0, zeros(d), fx)
  p0 = node0[:logf] - 0.5 * dot(node0[:r])
  logu = p0 + log(rand())
  xminus = xplus = node0[:x]
  rminus = rplus = node0[:r]
  gradminus = gradplus = node0[:grad]
  j = 0
  n = 1
  s = true
  node = Dict()
  while s
    pm = rand() > 0.5 ? 1 : -1
    if pm == -1
      node = buildtree(xminus, rminus, gradminus, logu, pm, j, eps, p0, fx)
      xminus = node[:xminus]
      rminus = node[:rminus]
      gradminus = node[:gradminus]
    else
      node = buildtree(xplus, rplus, gradplus, logu, pm, j, eps, p0, fx)
      xplus = node[:xplus]
      rplus = node[:rplus]
      gradplus = node[:gradplus]
    end
    if node[:s] && rand() < node[:n] / n
      v[:] = node[:x]
    end
    n += node[:n]
    s = node[:s] && nouturn(xminus, xplus, rminus, rplus)
    j += 1
  end
  v.tune.alpha = node[:alpha]
  v.tune.nalpha = node[:nalpha]
  v
end

function buildtree(x::Vector, r::Vector, grad::Vector, logu::Real, pm::Integer,
           j::Integer, eps::Real, p0::Real, fx::Function)
  if j == 0
    node = leapfrog(x, r, pm * eps, grad, fx)
    p = node[:logf] - 0.5 * dot(node[:r])
    node[:n] = logu <= p ? 1 : 0
    node[:s] = logu < p + 1000.0
    node[:xminus] = node[:xplus] = node[:x]
    node[:rminus] = node[:rplus] = node[:r]
    node[:gradminus] = node[:gradplus] = node[:grad]
    node[:alpha] = min(1.0, exp(p - p0))
    node[:nalpha] = 1
  else
    node = buildtree(x, r, grad, logu, pm, j-1, eps, p0, fx)
    if node[:s]
      if pm == -1
        node2 = buildtree(node[:xminus], node[:rminus], node[:gradminus],
                          logu, pm, j-1, eps, p0, fx)
        node[:xminus] = node2[:xminus]
        node[:rminus] = node2[:rminus]
        node[:gradminus] = node2[:gradminus]
      else
        node2 = buildtree(node[:xplus], node[:rplus], node[:gradplus],
                          logu, pm, j-1, eps, p0, fx)
        node[:xplus] = node2[:xplus]
        node[:rplus] = node2[:rplus]
        node[:gradplus] = node2[:gradplus]
      end
      if rand() < node2[:n] / (node[:n] + node2[:n])
        node[:x] = node2[:x]
      end
      node[:s] = node[:s] && node2[:s] &&
                 nouturn(node[:xminus], node[:xplus], node[:rminus], node[:rplus])
      node[:n] += node2[:n]
      node[:alpha] += node2[:alpha]
      node[:nalpha] += node2[:nalpha]
    end
  end
  node
end

function nouturn(xminus, xplus, rminus, rplus)
  val = xplus - xminus
  dot(val, rminus) >= 0 && dot(val, rplus) >= 0
end
