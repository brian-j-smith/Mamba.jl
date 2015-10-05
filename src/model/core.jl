#################### Core Model Functionality ####################

#################### Constructors ####################

function Model(; iter::Integer=0, burnin::Integer=0, chain::Integer=1,
               samplers::Vector{Sampler}=Sampler[], nodes...)
  nodedict = Dict{Symbol,Any}()
  for (key, value) in nodes
    isa(value, AbstractDependent) || error("nodes must be Dependent types")
    node = deepcopy(value)
    node.symbol = key
    nodedict[key] = node
  end
  m = Model(nodedict, Symbol[], Sampler[], Vector{Float64}[], iter, burnin,
            chain, false, false)
  g = graph(m)
  V = vertices(g)
  lookup = Dict{Symbol,Int}()
  for v in V
    setindex!(lookup, v.index, v.key)
  end
  m.dependents = intersect(tsort(g), keys(m, :dependent))
  for key in m.dependents
    targets = gettargets(V[lookup[key]], g, m)
    m[key].targets = intersect(m.dependents, targets)
  end
  setsamplers!(m, samplers)
end


#################### Indexing ####################

function Base.getindex(m::Model, key::Symbol)
  m.nodes[key]
end

function Base.setindex!(m::Model, values::Dict, nkeys::Vector{Symbol})
  for key in nkeys
    m[key][:] = values[key]
  end
end

function Base.setindex!(m::Model, value, nkeys::Vector{Symbol})
  length(nkeys) == 1 || throw(BoundsError())
  m[nkeys[1]][:] = value
end

function Base.setindex!(m::Model, value, nkey::Symbol)
  m[nkey][:] = value
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
      if isa(node, AbstractDependent) && length(node.monitor) > 0
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


#################### Auxiliary Functions ####################

function names(m::Model, monitoronly::Bool)
  values = AbstractString[]
  for key in keys(m, :dependent)
    node = m[key]
    lnames = unlist(node, names(node))
    v = monitoronly ? lnames[node.monitor] : vec(lnames)
    append!(values, v)
  end
  values
end

function names(m::Model, nkeys::Vector{Symbol})
  values = AbstractString[]
  for key in nkeys
    node = m[key]
    lnames = unlist(node, names(node))
    append!(values, vec(lnames))
  end
  values
end
