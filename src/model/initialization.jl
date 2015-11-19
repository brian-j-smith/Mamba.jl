#################### Model Initialization ####################

function reset!(m::Model, iter::Integer)
  m.iter = iter
  if iter == 0
    for s in m.samplers
      s.tune["sampler"] = nothing
    end
  end
  m
end

function setinits!(m::Model, inits::Dict{Symbol, Any})
  m.hasinputs || throw(ArgumentError("inputs must be set before inits"))
  reset!(m, 0)
  for key in keys(m, :dependent)
    node = m[key]
    if isa(node, AbstractStochastic)
      haskey(inits, key) ||
        throw(ArgumentError("missing initial value for node : $key"))
      setinits!(node, m, inits[key])
    else
      setinits!(node, m)
    end
  end
  m.hasinits = true
  m
end

function setinits!(m::Model, inits::Vector{Dict{Symbol, Any}})
  n = length(inits)
  m.states = Array(Vector{Float64}, n)
  for i in n:-1:1
    setinits!(m, inits[i])
    m.states[i] = unlist(m)
  end
  m
end

function setinputs!(m::Model, inputs::Dict{Symbol, Any})
  for key in keys(m, :input)
    haskey(inputs, key) ||
      throw(ArgumentError("missing inputs for node : $key"))
    isa(inputs[key], AbstractDependent) &&
      throw(ArgumentError("inputs cannot be Dependent types"))
    m.nodes[key] = deepcopy(inputs[key])
  end
  m.hasinputs = true
  m
end

function setsamplers!(m::Model, samplers::Vector{Sampler})
  m.samplers = deepcopy(samplers)
  for sampler in m.samplers
    sampler.targets = keys(m, :target, sampler.params)
  end
  m
end

function setstate!(m::Model, state::Vector{Float64}, iter::Integer)
  reset!(m, iter)
  relist!(m, state)
end
