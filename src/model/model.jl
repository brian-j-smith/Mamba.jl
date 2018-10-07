#################### Core Model Functionality ####################

#################### Constructors ####################

function Model(; iter::Integer=0, burnin::Integer=0,
               samplers::Vector{Sampler}=Sampler[], nodes...)
  nodedict = Dict{Symbol, Any}()
  for (key, value) in nodes
    isa(value, AbstractDependent) ||
      throw(ArgumentError("nodes are not all Dependent types"))
    node = deepcopy(value)
    node.symbol = key
    nodedict[key] = node
  end
  m = Model(nodedict, Sampler[], ModelState[], iter, burnin, false, false)
  dag = ModelGraph(m)
  dependentkeys = keys(m, :dependent)
  terminalkeys = keys(m, :stochastic)
  for v in vertices(dag.graph)
    vkey = dag.keys[v]
    if vkey in dependentkeys
      m[vkey].targets = intersect(dependentkeys,
                                  gettargets(dag, v, terminalkeys))
    end
  end
  setsamplers!(m, samplers)
end


#################### Indexing ####################

Base.getindex(m::Model, nodekey::Symbol) = m.nodes[nodekey]


function Base.setindex!(m::Model, value, nodekey::Symbol)
  node = m[nodekey]
  if isa(node, AbstractDependent)
    node.value = value
  else
    m.nodes[nodekey] = convert(typeof(node), value)
  end
end

function Base.setindex!(m::Model, values::Dict, nodekeys::Vector{Symbol})
  for key in nodekeys
    m[key] = values[key]
  end
end

function Base.setindex!(m::Model, value, nodekeys::Vector{Symbol})
  length(nodekeys) == 1 || throw(BoundsError())
  m[first(nodekeys)] = value
end


Base.keys(m::Model) = collect(keys(m.nodes))

function Base.keys(m::Model, ntype::Symbol, at...)
  ntype == :block       ? keys_block(m, at...) :
  ntype == :all         ? keys_all(m) :
  ntype == :assigned    ? keys_assigned(m) :
  ntype == :dependent   ? keys_dependent(m) :
  ntype == :independent ? keys_independent(m) :
  ntype == :input       ? keys_independent(m) :
  ntype == :logical     ? keys_logical(m) :
  ntype == :monitor     ? keys_monitor(m) :
  ntype == :output      ? keys_output(m) :
  ntype == :source      ? keys_source(m, at...) :
  ntype == :stochastic  ? keys_stochastic(m) :
  ntype == :target      ? keys_target(m, at...) :
    throw(ArgumentError("unsupported node type $ntype"))
end

function keys_all(m::Model)
  values = Symbol[]
  for key in keys(m)
    node = m[key]
    if isa(node, AbstractDependent)
      push!(values, key)
      append!(values, node.sources)
    end
  end
  unique(values)
end

function keys_assigned(m::Model)
  if m.hasinits
    values = keys(m)
  else
    values = Symbol[]
    for key in keys(m)
      if !isa(m[key], AbstractDependent)
        push!(values, key)
      end
    end
  end
  values
end

function keys_block(m::Model, block::Integer=0)
  block == 0 ? keys_block0(m) : m.samplers[block].params
end

function keys_block0(m::Model)
  values = Symbol[]
  for sampler in m.samplers
    append!(values, sampler.params)
  end
  unique(values)
end

function keys_dependent(m::Model)
  values = Symbol[]
  for key in keys(m)
    if isa(m[key], AbstractDependent)
      push!(values, key)
    end
  end
  intersect(tsort(m), values)
end

function keys_independent(m::Model)
  deps = Symbol[]
  for key in keys(m)
    if isa(m[key], AbstractDependent)
      push!(deps, key)
    end
  end
  setdiff(keys(m, :all), deps)
end

function keys_logical(m::Model)
  values = Symbol[]
  for key in keys(m)
    if isa(m[key], AbstractLogical)
      push!(values, key)
    end
  end
  values
end

function keys_monitor(m::Model)
  values = Symbol[]
  for key in keys(m)
    node = m[key]
    if isa(node, AbstractDependent) && !isempty(node.monitor)
      push!(values, key)
    end
  end
  values
end

function keys_output(m::Model)
  values = Symbol[]
  dag = ModelGraph(m)
  for v in vertices(dag.graph)
    vkey = dag.keys[v]
    if isa(m[vkey], AbstractStochastic) && !any_stochastic(dag, v, m)
      push!(values, vkey)
    end
  end
  values
end

keys_source(m::Model, nodekey::Symbol) = m[nodekey].sources

function keys_source(m::Model, nodekeys::Vector{Symbol})
  values = Symbol[]
  for key in nodekeys
    append!(values, m[key].sources)
  end
  unique(values)
end

function keys_stochastic(m::Model)
  values = Symbol[]
  for key in keys(m)
    if isa(m[key], AbstractStochastic)
      push!(values, key)
    end
  end
  values
end

function keys_target(m::Model, block::Integer=0)
  block == 0 ? keys_target0(m) : m.samplers[block].targets
end

function keys_target0(m::Model)
  values = Symbol[]
  for sampler in m.samplers
    append!(values, sampler.targets)
  end
  intersect(keys(m, :dependent), values)
end

keys_target(m::Model, nodekey::Symbol) = m[nodekey].targets

function keys_target(m::Model, nodekeys::Vector{Symbol})
  values = Symbol[]
  for key in nodekeys
    append!(values, m[key].targets)
  end
  intersect(keys(m, :dependent), values)
end


#################### Display ####################

function Base.show(io::IO, m::Model)
  showf(io, m, Base.show)
end

function showall(io::IO, m::Model)
  showf(io, m, Base.showall)
end

function showf(io::IO, m::Model, f::Function)
  print(io, "Object of type \"$(summary(m))\"\n")
  width = displaysize()[2] - 1
  for node in keys(m)
    print(io, string("-"^width, "\n", node, ":\n"))
    f(io, m[node])
    println(io)
  end
end


#################### Auxiliary Functions ####################

function names(m::Model, monitoronly::Bool)
  values = AbstractString[]
  for key in keys(m, :dependent)
    nodenames = names(m, key)
    v = monitoronly ? nodenames[m[key].monitor] : nodenames
    append!(values, v)
  end
  values
end

function names(m::Model, nodekey::Symbol)
  node = m[nodekey]
  unlist(node, names(node))
end

function names(m::Model, nodekeys::Vector{Symbol})
  values = AbstractString[]
  for key in nodekeys
    append!(values, names(m, key))
  end
  values
end
