#################### Missing Values Sampler ####################

#################### Types and Constructors ####################

type StochasticIndices
  dims::Tuple
  value::Vector{Int}
  distr::Vector{Int}
end


#################### Sampler Constructor ####################

function MISS(params::Vector{Symbol})
  Sampler(params,
    quote
      tunepar = tune(model, block)
      value = Dict{Symbol, Any}()
      initialize = tunepar["sampler"] == nothing
      if initialize
        tunepar["sampler"] = Dict{Symbol, StochasticIndices}()
      end
      for key in keys(model, :block, block)
        node = model[key]
        x = node[:]
        if initialize
          tunepar["sampler"][key] = findmissing(node)
        end
        inds = tunepar["sampler"][key]
        x[inds.value] = rand(node, inds)
        value[key] = x
      end
      value
    end,
    Dict{AbstractString, Any}("sampler" => nothing)
  )
end


#################### Sampling Functions ####################

function findmissing(s::AbstractStochastic)
  findmissing(s.distr, s.value)
end

function findmissing(d::Distribution, v)
  StochasticIndices((), find(isnan(v)), Int[])
end

function findmissing(D::Array{UnivariateDistribution}, v::Array)
  inds = find(isnan(v))
  StochasticIndices(dims(D), inds, inds)
end

function findmissing(D::Array{MultivariateDistribution}, v::Array)
  Mvalue = falses(v)
  Mdistr = falses(D)
  for sub in CartesianRange(size(D))
    n = length(D[sub])
    for i in 1:n
      if isnan(v[sub, i])
        Mvalue[sub, i] = Mdistr[sub] = true
      end
    end
  end
  StochasticIndices(dims(D), find(Mvalue), find(Mdistr))
end


rand(s::AbstractStochastic, inds::StochasticIndices) = rand_sub(s.distr, inds)

function rand_sub(d::Distribution, inds::StochasticIndices)
  if isempty(inds.value)
    Float64[]
  else
    x = rand(d)
    Float64[x[i] for i in inds.value]
  end
end

function rand_sub(D::Array{UnivariateDistribution}, inds::StochasticIndices)
  Float64[rand(d) for d in D[inds.value]]
end

function rand_sub(D::Array{MultivariateDistribution}, inds::StochasticIndices)
  X = Array(Float64, inds.dims)
  for i in inds.distr
    d = D[i]
    X[ind2sub(D, i)..., 1:length(d)] = rand(d)
  end
  X[inds.value]
end
