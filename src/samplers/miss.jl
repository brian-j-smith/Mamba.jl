#################### Missing Value Sampler ####################

function MISS(params::Vector{Symbol})
  Sampler(params,
    quote
      value = (Symbol => Any)[]
      sampler = model.samplers[block]
      for key in sampler.params
        node = model[key]
        v = deepcopy(node.value)
        if !sampler.tune["initialized"]
          sampler.tune["missing"][key] = find(isnan(node))
        end
        missing = sampler.tune["missing"][key]
        if isa(node.distr, Array)
          for i in missing
            v[i] = rand(node.distr[i])
          end
        elseif length(missing) > 0
          pred = rand(node.distr)
          for i in missing
            v[i] = pred[i]
          end
        end
        value[key] = v
      end
      sampler.tune["initialized"] = true
      value
    end,
    ["initialized" => false, "missing" => Dict{Symbol,Vector{Int}}()]
  )
end
