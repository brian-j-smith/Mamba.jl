#################### Core Model Functionality ####################

#################### Constructors ####################

function Model(; iter::Integer=0, burnin::Integer=0, chain::Integer=1,
               samplers::Vector{Sampler}=Sampler[], nodes...)
  nodedict = Dict{Symbol,Any}()
  for (key, value) in nodes
    isa(value, AbstractDependent) ||
      throw(ArgumentError("nodes are not all Dependent types"))
    node = deepcopy(value)
    node.symbol = key
    nodedict[key] = node
  end
  m = Model(nodedict, Symbol[], Sampler[], Vector{Float64}[], iter, burnin,
            chain, false, false)
  g = graph(m)
  m.dependents = intersect(tsort(g), keys(m, :dependent))
  for v in vertices(g)
    if v.key in m.dependents
      m[v.key].targets = intersect(m.dependents, gettargets(v, g, m))
    end
  end
  setsamplers!(m, samplers)
end


#################### Indexing ####################

function Base.getindex(m::Model, key::Symbol)
  m.nodes[key]
end

function Base.setindex!(m::Model, values::Dict, nodekeys::Vector{Symbol})
  for key in nodekeys
    m[key][:] = values[key]
  end
end

function Base.setindex!(m::Model, value, nodekeys::Vector{Symbol})
  length(nodekeys) == 1 || throw(BoundsError())
  m[nodekeys[1]][:] = value
end

function Base.setindex!(m::Model, value, nodekey::Symbol)
  m[nodekey][:] = value
end


#################### Base Methods ####################

function Base.keys(m::Model, ntype::Symbol=:assigned, block::Integer=0)
  values = Symbol[]
  if ntype == :all
    for key in keys(m.nodes)
      if isa(m[key], AbstractDependent)
        values = [values; key; m[key].sources]
      end
    end
    values = unique(values)
  elseif ntype == :assigned
    values = collect(keys(m.nodes))
  elseif ntype == :block
    blocks = block > 0 ? block : 1:length(m.samplers)
    for b in blocks
      append!(values, m.samplers[b].params)
    end
    values = unique(values)
  elseif ntype == :dependent
    for key in keys(m.nodes)
      if isa(m[key], AbstractDependent)
        push!(values, key)
      end
    end
  elseif ntype == :independent || ntype == :input
    values = setdiff(keys(m, :all), keys(m, :dependent))
  elseif ntype == :logical
    for key in keys(m.nodes)
      if isa(m[key], AbstractLogical)
        push!(values, key)
      end
    end
  elseif ntype == :monitor
    for key in keys(m.nodes)
      node = m[key]
      if isa(node, AbstractDependent) && !isempty(node.monitor)
        push!(values, key)
      end
    end
  elseif ntype == :output
    g = graph(m)
    for v in vertices(g)
      if isa(m[v.key], AbstractStochastic) && !any_stochastic(v, g, m)
        push!(values, v.key)
      end
    end
  elseif ntype == :stochastic
    for key in keys(m.nodes)
      if isa(m[key], AbstractStochastic)
        push!(values, key)
      end
    end
  else
    error("unsupported node type $ntype")
  end
  values
end


#################### Display ####################

function Base.show(io::IO, m::Model)
  showf(io, m, Base.show)
end

function Base.showall(io::IO, m::Model)
  showf(io, m, Base.showall)
end

function showf(io::IO, m::Model, f::Function)
  print(io, "Object of type \"$(summary(m))\"\n")
  width = Base.tty_size()[2] - 1
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
