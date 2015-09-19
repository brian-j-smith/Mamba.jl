#################### Missing Values Sampler ####################

#################### Types ####################

type StochasticIndices
  dims::Tuple
  value::Vector{Int}
  distr::Vector{Int}
end

function StochasticIndices(s::AbstractStochastic, valueinds::Vector{Int})
  StochasticIndices(size(s), valueinds, mapindices(s.distr, valueinds))
end


#################### Sampler Constructor ####################

function MISS(params::Vector{Symbol})
  Sampler(params,
    quote
      tunepar = tune(model, block)
      value = Dict{Symbol,Any}()
      initialize = tunepar["sampler"] == nothing
      if initialize
        tunepar["sampler"] = Dict{Symbol,StochasticIndices}()
      end
      for key in keys(model, :block, block)
        node = model[key]
        v = node[:]
        if initialize
          tunepar["sampler"][key] = StochasticIndices(node, find(isnan(node)))
        end
        inds = tunepar["sampler"][key]
        v[inds.value] = sample(node, inds)
        value[key] = v
      end
      value
    end,
    Dict{String,Any}("sampler" => nothing)
  )
end


#################### Sampling Functions ####################

mapindices(d::Distribution, valueinds::Vector{Int}) = valueinds

mapindices(D::Array{UnivariateDistribution}, valueinds::Vector{Int}) = valueinds

function mapindices(D::Array{MultivariateDistribution}, valueinds::Vector{Int})
  M = falses(D)
  valuedims = dim(D)
  for i in valueinds
    sub = ind2sub(valuedims, i)[1:ndims(M)]
    M[sub...] = true
  end
  find(M)
end

sample(s::AbstractStochastic, inds::StochasticIndices) = sample(s.distr, inds)

function sample(d::Distribution, inds::StochasticIndices)
  length(inds.value) > 0 ? rand(d)[inds.value] : Float64[]
end

function sample(D::Array{UnivariateDistribution}, inds::StochasticIndices)
  Float64[rand(d) for d in D[inds.distr]]
end

function sample(D::Array{MultivariateDistribution}, inds::StochasticIndices)
  X = Array(Float64, inds.dims)
  for i in inds.distr
    sub = ind2sub(size(D), i)
    X[sub..., :] = rand(D[sub...])
  end
  X[inds.value]
end
