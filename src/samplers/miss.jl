#################### Missing Value Sampler ####################

function MISS(params::Vector{Symbol})
  Sampler(params,
    quote
      value = (Symbol => Any)[]
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
      value
    end,
    ["sampler" => nothing]
  )
end
