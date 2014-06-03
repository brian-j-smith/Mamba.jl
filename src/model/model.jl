#################### MCMCModel Constructor ####################

function MCMCModel(; iter::Integer=0, burnin::Integer=0, chain::Integer=1,
           samplers::Vector{MCMCSampler}=MCMCSampler[], nodes...)
  nodedict = (Symbol => Any)[]
  for (key, value) in nodes
    isa(value, MCMCDependent) || error("nodes must be MCMCDependent types")
    node = deepcopy(value)
    node.symbol = key
    nodedict[key] = node
  end
  m = MCMCModel(nodedict, Symbol[], MCMCSampler[], iter, burnin, chain, false,
                false)
  g = graph(m)
  V = vertices(g)
  lookup = (Symbol => Integer)[]
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


#################### MCMCModel Base Methods ####################

function Base.getindex(m::MCMCModel, key::Symbol)
  m.nodes[key]
end

function Base.keys(m::MCMCModel, ntype::Symbol=:assigned, block::Integer=0)
  values = Symbol[]
  if ntype == :all
    for key in keys(m.nodes)
      if isa(m[key], MCMCDependent)
        values = [values, key, m[key].sources]
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
      if isa(m[key], MCMCDependent)
        push!(values, key)
      end
    end
  elseif ntype == :independent || ntype == :input
    values = setdiff(keys(m, :all), keys(m, :dependent))
  elseif ntype == :logical
    for key in keys(m.nodes)
      if isa(m[key], MCMCLogical)
        push!(values, key)
      end
    end
  elseif ntype == :monitor
    for key in keys(m.nodes)
      node = m[key]
      if isa(node, MCMCDependent) && any(node.monitor)
        push!(values, key)
      end
    end
  elseif ntype == :output
    g = graph(m)
    for v in vertices(g)
      if isa(m[v.key], MCMCStochastic) && !any_stochastic(v, g, m)
        push!(values, v.key)
      end
    end
  elseif ntype == :stochastic
    for key in keys(m.nodes)
      if isa(m[key], MCMCStochastic)
        push!(values, key)
      end
    end
  else
    error("unsupported node type $ntype")
  end
  values
end

function Base.setindex!(m::MCMCModel, values::Dict, nkeys::Vector{Symbol})
  for key in nkeys
    m[key][:] = values[key]
  end
end

function Base.setindex!(m::MCMCModel, value, nkeys::Vector{Symbol})
  length(nkeys) == 1 || throw(BoundsError())
  m[nkeys[1]][:] = value
end

function Base.setindex!(m::MCMCModel, value, nkey::Symbol)
  m[nkey][:] = value
end

function Base.show(io::IO, m::MCMCModel)
  showf(io, m, Base.show)
end

function Base.showall(io::IO, m::MCMCModel)
  showf(io, m, Base.showall)
end

function showf(io::IO, m::MCMCModel, f::Function)
  print(io, "Object of type \"$(summary(m))\"\n")
  width = Base.tty_cols() - 1
  for node in keys(m)
    print(io, string("-"^width, "\n", node, ":\n"))
    f(io, m[node])
    println(io)
  end
end


#################### MCMCModel Initialization Methods ####################

function names(m::MCMCModel, monitoronly::Bool)
  values = String[]
  for key in keys(m, :dependent)
    node = m[key]
    append!(values, names(node)[!monitoronly | node.monitor])
  end
  values
end

function names(m::MCMCModel, nkeys::Vector{Symbol})
  values = String[]
  for key in nkeys
    append!(values, names(m[key]))
  end
  values
end

function setinits!(m::MCMCModel, inits::Dict{Symbol,Any})
  m.iter = 0
  for key in m.dependents
    node = m[key]
    if isa(node, MCMCStochastic)
      setinits!(node, m, inits[key])
    else
      setinits!(node, m)
    end
  end
  m
end

function setinputs!(m::MCMCModel, inputs::Dict{Symbol,Any})
  for key in keys(m, :input)
    isa(inputs[key], MCMCDependent) &&
      error("inputs must not be MCMCDependent types")
    m.nodes[key] = deepcopy(inputs[key])
  end
  m.hasinputs = true
  m
end

function setsamplers!(m::MCMCModel, samplers::Vector{MCMCSampler})
  m.samplers = deepcopy(samplers)
  for i in 1:length(m.samplers)
    sampler = m.samplers[i]
    targets = mapreduce(key -> m[key].targets, vcat, sampler.params)
    sampler.targets = intersect(m.dependents, targets)
  end
  m
end

function settune!(m::MCMCModel, tune::Vector)
  for b in 1:length(m.samplers)
    m.samplers[b].tune = deepcopy(tune[b])
  end
end

function tune(m::MCMCModel, block::Integer=0)
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


#################### MCMCModel Simulation Methods ####################

function gradient(m::MCMCModel, block::Integer=0, transform::Bool=false;
           dtype::Symbol=:forward)
  x0 = unlist(m, block, transform)
  value = gradient!(m, x0, block, transform, dtype=dtype)
  relist!(m, x0, block, transform)
  value
end

function gradient{T<:Real}(m::MCMCModel, x::Vector{T}, block::Integer=0,
           transform::Bool=false; dtype::Symbol=:forward)
  x0 = unlist(m, block)
  value = gradient!(m, x, block, transform, dtype=dtype)
  relist!(m, x0, block)
  value
end

function gradient!{T<:Real}(m::MCMCModel, x::Vector{T}, block::Integer=0,
           transform::Bool=false; dtype::Symbol=:forward)
  f = x -> logpdf!(m, x, block, transform)
  gradient(f, x, dtype)
end

function logpdf(m::MCMCModel, block::Integer=0, transform::Bool=false)
  value = 0
  if block > 0
    sampler = m.samplers[block]
    params = sampler.params
    nkeys = [setdiff(params, sampler.targets), sampler.targets]
  else
    params = keys(m, :block)
    nkeys = m.dependents
  end
  for key in nkeys
    value += logpdf(m[key], transform && in(key, params))
    isfinite(value) || break
  end
  value
end

function logpdf{T<:Real}(m::MCMCModel, x::Vector{T}, block::Integer=0,
           transform::Bool=false)
  x0 = unlist(m, block)
  value = logpdf!(m, x, block, transform)
  relist!(m, x0, block)
  value
end

function logpdf!{T<:Real}(m::MCMCModel, x::Vector{T}, block::Integer=0,
           transform::Bool=false)
  value = 0
  if block > 0
    sampler = m.samplers[block]
    params = sampler.params
    targets = sampler.targets
  else
    params = keys(m, :block)
    targets = m.dependents
  end
  m[params] = relist(m, x, params, transform)
  for key in setdiff(params, targets)
    value += logpdf(m[key], transform)
    isfinite(value) || return value
  end
  for key in targets
    update!(m[key], m)
    value += logpdf(m[key], transform && in(key, params))
    isfinite(value) || return value
  end
  value
end

function relist{T<:Real}(m::MCMCModel, values::Vector{T}, block::Integer=0,
           transform::Bool=false)
  relist(m, values, keys(m, :block, block), transform)
end

function relist{T<:Real}(m::MCMCModel, values::Vector{T}, nkeys::Vector{Symbol},
           transform::Bool=false)
  f =  transform ? invlink : identity
  x = (Symbol => Any)[]
  j = 0
  for key in nkeys
    node = m[key]
    n = length(node)
    x[key] = f(node, values[j+(1:n)])
    j += n
  end
  j == length(values) || throw(ErrorException("argument dimensions must match"))
  x
end

function relist!{T<:Real}(m::MCMCModel, values::Vector{T}, block::Integer=0,
           transform::Bool=false)
  nkeys = keys(m, :block, block)
  m[nkeys] = relist(m, values, nkeys, transform)
  update!(m, block)
end

function relist!{T<:Real}(m::MCMCModel, values::Vector{T},
           nkeys::Vector{Symbol}, transform::Bool=false)
  m[nkeys] = relist(m, values, nkeys, transform)
  update!(m)
end

function simulate!(m::MCMCModel, block::Integer=0)
  if block > 0
    blocks = block
  else
    m.iter += 1
    blocks = 1:length(m.samplers)
  end
  for b in blocks
    sampler = m.samplers[b]
    value = sampler.eval(m, b)
    if value != nothing
      m[sampler.params] = value
      update!(m, b)
    end
  end
  m
end

function unlist(m::MCMCModel, block::Integer=0, transform::Bool=false)
  unlist(m, keys(m, :block, block), transform)
end

function unlist(m::MCMCModel, monitoronly::Bool)
  values = VariateType[]
  for key in keys(m, :dependent)
    node = m[key]
    idx = find(!monitoronly | node.monitor)
    if length(idx) > 0
      append!(values, node[idx])
    end
  end
  values
end

function unlist(m::MCMCModel, nkeys::Vector{Symbol}, transform::Bool=false)
  f = transform ? link : identity
  N = map(key -> length(m[key]), nkeys)
  values = Array(VariateType, sum(N))
  i = 0
  for k in 1:length(nkeys)
    node = m[nkeys[k]]
    n = N[k]
    values[i+(1:n)] = f(node, node.value)
    i += n
  end
  values
end

function update!(m::MCMCModel, block::Integer=0)
  nkeys = block > 0 ? m.samplers[block].targets : m.dependents
  for key in nkeys
    update!(m[key], m)
  end
  m
end
