#################### Discrete Gibbs Sampler ####################

#################### Types and Constructors ####################

typealias DGSUnivariateDistribution
          Union{Bernoulli, Binomial, Categorical, DiscreteUniform,
                Hypergeometric, NoncentralHypergeometric}

type DGSTune
  support::Matrix{Real}
  probs::Vector{Float64}
end

type DGSVariate <: VectorVariate
  value::Vector{Float64}
  tune::DGSTune

  DGSVariate(x::Vector{Float64}, tune::DGSTune) = new(x, tune)
end

function DGSVariate(x::Vector{Float64}, tune=nothing)
  tune = DGSTune(
    Array(Float64, 0, 0),
    Array(Float64, 0)
  )
  DGSVariate(x, tune)
end


#################### Sampler Constructor ####################

function DGS(params::Vector{Symbol})
  Sampler(params,
    quote
      tunepar = tune(model, block)
      x = unlist(model, block)
      offset = 0
      for key in keys(model, :block, block)

        sim = function(inds, d, logf)
          v = DGSVariate(x[offset + inds], tunepar["sampler"])
          dgs!(v, d, logf)
          x[offset + inds] = v
          tunepar["sampler"] = v.tune
        end

        logf = function(inds, value)
          x[offset + inds] = value
          logpdf!(model, x, block)
        end

        node = model[key]
        DGS_sub!(node.distr, sim, logf)
        offset += length(node)
      end
      relist(model, x, block)
    end,
    Dict{AbstractString,Any}("sampler" => nothing)
  )
end

function DGS_sub!(d::UnivariateDistribution, sim::Function, logf::Function)
  inds = [1]
  sim(inds, d, x -> logf(inds, x))
end

function DGS_sub!(D::Array{UnivariateDistribution}, sim::Function,
                  logf::Function)
  for i in 1:length(D)
    inds = [i]
    sim(inds, D[i], x -> logf(inds, x))
  end
end

function DGS_sub!(d, sim::Function, logf::Function)
  throw(ArgumentError("unsupported distribution structure $(typeof(d))"))
end


#################### Sampling Functions ####################

function dgs!{T<:Real}(v::DGSVariate, support::Matrix{T}, logf::Function)
  n = size(support, 1)
  probs = Array(Float64, n)
  psum = 0.0
  for i in 1:n
    x = vec(support[i,:])
    value = exp(logf(x))
    probs[i] = value
    psum += value
  end
  if psum > 0
    probs /= psum
  else
    probs[:] = 1 / n
  end
  v[:] = support[rand(Categorical(probs)), :]
  v.tune.support = support
  v.tune.probs = probs
  v
end

function dgs!{T<:Real}(v::DGSVariate, support::Matrix{T},
                       probs::Vector{Float64})
  size(support, 1) == length(probs) ||
    throw(ArgumentError("number of support rows and probs length differ"))

  v[:] = support[rand(Categorical(probs)), :]
  v.tune.support = support
  v.tune.probs = probs
  v
end

function dgs!(v::DGSVariate, d::DGSUnivariateDistribution, logf::Function)
  S = support(d)
  dgs!(v, reshape(S, length(S), 1), logf)
end

function dgs!(v::DGSVariate, d::Distribution, logf::Function)
  throw(ArgumentError("unsupported distribution $(typeof(d))"))
end
