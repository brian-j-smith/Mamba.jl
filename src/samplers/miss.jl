#################### Missing Values Sampler ####################

#################### Sampler Constructor ####################

function MISS(params::Vector{Symbol})
  Sampler(params,
    quote
      value = Dict{Symbol,Any}()
      tunepar = tune(model, block)
      initialize = tunepar["sampler"] == nothing
      if initialize
        tunepar["sampler"] = Dict{Symbol,Vector{Int}}()
      end
      for key in keys(model, :block, block)
        node = model[key]
        v = copy(node.value)
        if initialize
          tunepar["sampler"][key] = find(isnan(node))
        end
        missing = tunepar["sampler"][key]
        v[missing] = sample(node.distr, missing)
        value[key] = v
      end
      value
    end,
    Dict{String,Any}("sampler" => nothing)
  )
end


#################### Sampling Functions ####################

function sample(d::Distribution, inds::Vector{Int})
  length(inds) > 0 ? rand(d)[inds] : Float64[]
end

function sample(D::Array{UnivariateDistribution}, inds::Vector{Int})
  Float64[rand(d) for d in D[inds]]
end

function sample(D::Array{MultivariateDistribution}, inds::Vector{Int})
  X = Array(Float64, dim(D))
  ID = falses(D)
  for i in inds
    sub = ind2sub(X, i)[1:ndims(D)]
    ID[sub...] = true
  end
  for sub in CartesianRange(size(D))
    if ID[sub]
      X[sub, :] = rand(D[sub])
    end
  end
  X[inds]
end
