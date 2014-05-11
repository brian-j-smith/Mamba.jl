#################### Missing Value Sampler ####################

function MISS{T<:String}(params::Vector{T})
  length(params) == 1 || error("must specify a single node to sample")
  MCMCSampler(params,
    quote
      sampler = model.samplers[block]
      node = model[sampler.params[1]]
      value = deepcopy(node.value)
      if !sampler.tune["initialized"]
        sampler.tune["missing"] = find(isnan(node))
      end
      missing = sampler.tune["missing"]
      if isa(node.distr, Array)
        for i in missing
          value[i] = rand(node.distr[i])
        end
      elseif length(missing) > 0
        pred = rand(node.distr)
        for i in missing
          value[i] = pred[i]
        end
      end
      value
    end,
    ["initialized" => false, "missing" => Int[]]
  )
end
