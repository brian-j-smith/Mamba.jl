#################### Missing Values Sampler ####################

#################### Types and Constructors ####################

type MISSTune
  dims::Tuple
  valueinds::Vector{Int}
  distrinds::Vector{Int}
end

function MISSTune(s::AbstractStochastic)
  MISSTune(s.distr, s.value)
end

function MISSTune(d::Distribution, v)
  MISSTune((), find(isnan(v)), Int[])
end

function MISSTune(D::Array{UnivariateDistribution}, v::Array)
  inds = find(isnan(v))
  MISSTune(dims(D), inds, inds)
end

function MISSTune(D::Array{MultivariateDistribution}, v::Array)
  isvalueinds = falses(v)
  isdistrinds = falses(D)
  for sub in CartesianRange(size(D))
    n = length(D[sub])
    for i in 1:n
      if isnan(v[sub, i])
        isvalueinds[sub, i] = isdistrinds[sub] = true
      end
    end
  end
  MISSTune(dims(D), find(isvalueinds), find(isdistrinds))
end


#################### Sampler Constructor ####################

function MISS(params::ElementOrVector{Symbol})
  params = asvec(params)
  samplerfx = function(model::Model, block::Integer)
    tunepar = tune(model, block)
    value = Dict{Symbol, Any}()
    isfirstiter = model.iter == 1
    for key in params
      node = model[key]
      x = node[:]
      if isfirstiter
        tunepar[key] = MISSTune(node)
      end
      miss = tunepar[key]
      x[miss.valueinds] = rand(node, miss)
      value[key] = x
    end
    value
  end
  Sampler(params, samplerfx, Dict{Symbol, MISSTune}())
end


#################### Sampling Functions ####################

rand(s::AbstractStochastic, miss::MISSTune) = rand_sub(s.distr, miss)

function rand_sub(d::Distribution, miss::MISSTune)
  if isempty(miss.valueinds)
    Float64[]
  else
    x = rand(d)
    Float64[x[i] for i in miss.valueinds]
  end
end

function rand_sub(D::Array{UnivariateDistribution}, miss::MISSTune)
  Float64[rand(d) for d in D[miss.valueinds]]
end

function rand_sub(D::Array{MultivariateDistribution}, miss::MISSTune)
  X = Array(Float64, miss.dims)
  for i in miss.distrinds
    d = D[i]
    X[ind2sub(D, i)..., 1:length(d)] = rand(d)
  end
  X[miss.valueinds]
end
