#################### Discrete Gibbs Sampler ####################

#################### Types ####################

typealias DGSUnivariateDistribution
          Union{Bernoulli, Binomial, Categorical, DiscreteUniform,
                Hypergeometric, NoncentralHypergeometric}

type DGSTune
  support::Vector
  probs::Vector{Float64}
end

type DGSVariate <: VectorVariate
  value::Vector{Float64}
  tune::DGSTune

  DGSVariate(x::Vector{Float64}, tune::DGSTune) = new(x, tune)
end

function DGSVariate(x::Vector{Float64}, tune=nothing)
  tune = DGSTune(
    Array(Any, 0),
    Array(Float64, 0)
  )
  DGSVariate(x, tune)
end


#################### Sampler Constructor ####################

function DGS(params::Vector{Symbol})
  Sampler(params,
    quote
      x = unlist(model, block)
      offset = 0
      for key in keys(model, :block, block)

        sim = function(inds, d, logf)
          x[offset + inds], _ = dgs(d, logf)
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
  sim(1, d, x -> logf(1, x))
end

function DGS_sub!(D::Array{UnivariateDistribution}, sim::Function,
                  logf::Function)
  for i in 1:length(D)
    sim(i, D[i], x -> logf(i, x))
  end
end

function DGS_sub!(d, sim::Function, logf::Function)
  throw(ArgumentError("unsupported distribution structure $(typeof(d))"))
end


#################### Sampling Functions ####################

function dgs!(v::DGSVariate, support::Vector, logf::Function)
  v[:], probs = dgs(support, logf)
  v.tune.support = support
  v.tune.probs = probs
  v
end

function dgs!(v::DGSVariate, support::Vector, probs::Vector{Float64})
  length(support) == length(probs) ||
    throw(ArgumentError("lengths of support and probs differ"))

  v[:] = support[rand(Categorical(probs))]
  v.tune.support = support
  v.tune.probs = probs
  v
end

function dgs(support::AbstractVector, logf::Function)
  n = length(support)
  p = Array(Float64, n)
  psum = 0.0
  for i in 1:n
    value = exp(logf(support[i]))
    p[i] = value
    psum += value
  end
  if psum > 0
    p /= psum
  else
    p[:] = 1 / n
  end
  support[rand(Categorical(p))], p
end

function dgs(d::DGSUnivariateDistribution, logf::Function)
  dgs(support(d), logf)
end

function dgs(d::Distribution, logf::Function)
  throw(ArgumentError("unsupported distribution $(typeof(d))"))
end
