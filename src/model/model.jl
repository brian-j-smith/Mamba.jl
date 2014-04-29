#################### MCMCModel Constructor ####################

function MCMCModel(; iter::Integer=0, burnin::Integer=0, chain::Integer=1,
           samplers::Vector{MCMCSampler}=MCMCSampler[], nodes...)
  nodedict = Dict{String,Any}()
  for (arg, value) in nodes
    isa(value, MCMCDependent) || error("nodes must be MCMCDependent types")
    node = deepcopy(value)
    key = string(arg)
    node.names = names(node, key)
    nodedict[key] = node
  end
  m = MCMCModel(nodedict, String[], MCMCSampler[], iter, burnin, chain, false,
                false)
  setsamplers!(m, samplers)
end


#################### MCMCModel Base Methods ####################

function Base.getindex(m::MCMCModel, key::String)
  m.nodes[key]
end

function Base.keys(m::MCMCModel, ntype::Symbol=:assigned, block::Integer=0)
  values = String[]
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
  elseif ntype == :dep
    for key in keys(m.nodes)
      if isa(m[key], MCMCDependent)
        push!(values, key)
      end
    end
  elseif ntype == :indep || ntype == :input
    values = setdiff(keys(m, :all), keys(m, :dep))
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
  elseif ntype == :terminal
    g = graph(m)
    for v in vertices(g)
      if isa(m[v.label], MCMCStochastic) && !any_stochastic(v, g, m)
        push!(values, v.label)
      end
    end
  elseif ntype == :stochastic
    for key in keys(m.nodes)
      if isa(m[key], MCMCStochastic)
        push!(values, key)
      end
    end
  else
    error("unsupported ntype")
  end
  values
end

function Base.setindex!{T<:String}(m::MCMCModel, values::Dict, nkeys::Vector{T})
  for key in nkeys
    m[key][:] = values[key]
  end
end

function Base.setindex!{T<:String}(m::MCMCModel, value, nkeys::Vector{T})
  length(nkeys) == 1 || throw(BoundsError())
  m[nkeys[1]][:] = value
end

function Base.setindex!{T<:String}(m::MCMCModel, value, key::T)
  m[key][:] = value
end

function Base.show(io::IO, m::MCMCModel)
  print(io, "Object of type \"$(summary(m))\"\n")
  for node in keys(m)
    print(io, string("-"^79, "\n", node, ":\n"))
    show(io, m[node])
    print(io, "\n")
  end
end

function Base.showall(io::IO, m::MCMCModel)
  print(io, "Object of type \"$(summary(m))\"\n")
  for node in keys(m)
    print(io, string("-"^79, "\n", node, ":\n"))
    showall(io, m[node])
    print(io, "\n")
  end
end


#################### MCMCModel Initialization Methods ####################

function names(m::MCMCModel, monitoronly::Bool)
  values = String[]
  for key in keys(m, :dep)
    node = m[key]
    append!(values, node.names[!monitoronly | node.monitor])
  end
  values
end

function names{T<:String}(m::MCMCModel, nkeys::Vector{T})
  values = String[]
  for key in nkeys
    append!(values, m[key].names)
  end
  values
end

function setinits!{T<:String}(m::MCMCModel, inits::Dict{T,Any})
  for key in m.targets
    node = m[key]
    if isa(node, MCMCStochastic)
      setinits!(m[key], m, inits[key])
    else
      setinits!(m[key], m)
    end
  end
  m
end

function setinputs!{T<:String}(m::MCMCModel, inputs::Dict{T,Any})
  for key in keys(m, :input)
    isa(inputs[key], MCMCDependent) &&
      error("inputs must not be MCMCDependent types")
    m.nodes[key] = deepcopy(inputs[key])
  end
  m.hasinputs = true
  m
end

function setsamplers!(m::MCMCModel, samplers::Vector{MCMCSampler})
  g = graph(m)
  m.targets = intersect(tsort(g), keys(m, :dep))
  m.samplers = deepcopy(samplers)
  V = vertices(g)
  lookup = Dict{String,Integer}()
  for v in V
    setindex!(lookup, v.index, v.label)
  end
  for i in 1:length(m.samplers)
    sampler = m.samplers[i]
    targets = String[]
    for key in sampler.params
      append!(targets, gettargets(V[lookup[key]], g, m))
    end
    sampler.targets = intersect(m.targets, targets)
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

function gradient(m::MCMCModel, block::Integer=0, transform::Bool=false,
           dtype::Symbol=:central)
  x0 = unlist(m, block, transform)
  value = gradient!(m, x0, block, transform, dtype)
  relist!(m, x0, block, transform)
  value
end

function gradient{T<:Real}(m::MCMCModel, x::Vector{T}, block::Integer=0,
           transform::Bool=false, dtype::Symbol=:central)
  x0 = unlist(m, block)
  value = gradient!(m, x, block, transform, dtype)
  relist!(m, x0, block)
  value
end

function gradient!{T<:Real}(m::MCMCModel, x::Vector{T}, block::Integer=0,
           transform::Bool=false, dtype::Symbol=:central)
  f = x -> logpdf!(m, x, block, transform)
  gradient(f, x, dtype)
end

function logpdf(m::MCMCModel, block::Integer=0, transform::Bool=false)
  blocks = block > 0 ? block : 1:length(m.samplers)
  nkeys = String[]
  for b in blocks
    append!(nkeys, m.samplers[b].params)
    append!(nkeys, m.samplers[b].targets)
  end
  mapreduce(key -> logpdf(m[key], transform), +, unique(nkeys))
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
  nkeys = keys(m, :block, block)
  m[nkeys] = relist(m, x, nkeys, transform)
  if all(map(key -> insupport(m[key]), nkeys))
    update!(m, block)
    logpdf(m, block, transform)
  else
    -Inf
  end
end

function relist{T<:Real}(m::MCMCModel, values::Vector{T}, block::Integer=0,
           transform::Bool=false)
  relist(m, values, keys(m, :block, block), transform)
end

function relist{T<:Real,U<:String}(m::MCMCModel, values::Vector{T},
           nkeys::Vector{U}, transform::Bool=false)
  f =  transform ? invlink : identity
  x = Dict{String,Any}()
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

function relist!{T<:Real,U<:String}(m::MCMCModel, values::Vector{T},
           nkeys::Vector{U}, transform::Bool=false)
  m[nkeys] = relist(m, values, nkeys, transform)
  update!(m)
end

function simulate!(m::MCMCModel, block::Integer=0)
  blocks = block > 0 ? block : 1:length(m.samplers)
  for b in blocks
    sampler = m.samplers[b]
    m[sampler.params] = sampler.eval(m, b)
    update!(m, b)
  end
  m
end

function unlist(m::MCMCModel, block::Integer=0, transform::Bool=false)
  unlist(m, keys(m, :block, block), transform)
end

function unlist(m::MCMCModel, monitoronly::Bool)
  values = VariateType[]
  for key in keys(m, :dep)
    node = m[key]
    idx = find(!monitoronly | node.monitor)
    if length(idx) > 0
      append!(values, node[idx])
    end
  end
  values
end

function unlist{T<:String}(m::MCMCModel, nkeys::Vector{T}, transform::Bool=false)
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
  targets = block > 0 ? m.samplers[block].targets : m.targets
  for key in targets
    update!(m[key], m)
  end
  m
end
