#################### Model Initialization Methods ####################

function reset!(m::Model)
  m.iter = 0
  for s in m.samplers
    s.tune["sampler"] = nothing
  end
  m
end

function setinits!(m::Model, inits::Dict{Symbol,Any})
  reset!(m)
  for key in m.dependents
    node = m[key]
    if isa(node, Stochastic)
      haskey(inits, key) || error(string("missing inits for node :", key))
      setinits!(node, m, inits[key])
    else
      setinits!(node, m)
    end
  end
  m
end

function setinits!{T<:Real}(m::Model, inits::Vector{T})
  reset!(m)
  relist!(m, inits)
end

function setinputs!(m::Model, inputs::Dict{Symbol,Any})
  for key in keys(m, :input)
    haskey(inputs, key) || error(string("missing inputs for node :", key))
    isa(inputs[key], Dependent) && error("inputs cannot be Dependent types")
    m.nodes[key] = deepcopy(inputs[key])
  end
  m.hasinputs = true
  m
end

function setsamplers!(m::Model, samplers::Vector{Sampler})
  m.samplers = deepcopy(samplers)
  for i in 1:length(m.samplers)
    sampler = m.samplers[i]
    targets = mapreduce(key -> m[key].targets, vcat, sampler.params)
    sampler.targets = intersect(m.dependents, targets)
  end
  m
end
