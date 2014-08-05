#################### Model Initialization Methods ####################

function names(m::Model, monitoronly::Bool)
  values = String[]
  for key in keys(m, :dependent)
    node = m[key]
    lnames = link(node, names(node), false)
    v = monitoronly ? lnames[node.monitor] : vec(lnames)
    append!(values, v)
  end
  values
end

function names(m::Model, nkeys::Vector{Symbol})
  values = String[]
  for key in nkeys
    node = m[key]
    lnames = link(node, names(node), false)
    append!(values, vec(lnames))
  end
  values
end

function setinits!(m::Model, inits::Dict{Symbol,Any})
  m.iter = 0
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
  m.iter = 0
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

function settune!(m::Model, tune::Vector)
  for b in 1:length(m.samplers)
    m.samplers[b].tune = deepcopy(tune[b])
  end
end

function tune(m::Model, block::Integer=0)
  if block > 0
    values = m.samplers[block].tune
  else
    n = length(m.samplers)
    values = Array(Any, n)
    for i in 1:n
      values[i] = m.samplers[i].tune
    end
  end
  values
end
