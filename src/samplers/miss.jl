#################### Missing Value Sampler ####################

function MISS{T<:String}(params::Vector{T})
  length(params) == 1 || error("must specify a single node to sample")
  MCMCSampler(params,
    quote
      sampler = model.samplers[block]
      node = model[sampler.params[1]]
      if model.iter == 1
        sampler.tune["missing"] = isnan(node)
      end
      missing = sampler.tune["missing"]
      if isa(node.distr, Array)
        for i in find(missing)
          node[i] = rand(node.distr[i])
        end
      else
        node[missing] = rand(node.distr)[missing]
      end
      node.value
    end,
    (String => Any)["missing" => nothing]
  )
end
