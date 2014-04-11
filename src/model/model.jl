#################### MCMCModel Constructor ####################

function MCMCModel(; iter::Integer=0, burnin::Integer=0, chain::Integer=1,
           samplers::Vector{MCMCSampler}=MCMCSampler[], nodes...)
  nodedict = Dict{String,Any}()
  for (key, value) in nodes
    isa(value, MCMCDepNode) || error("nodes must be MCMCDepNode types")
    nodedict[string(key)] = deepcopy(value)
  end
  m = MCMCModel(nodedict, String[], samplers, iter, burnin, chain, false, false)
  g = graph(m)
  m.links = intersect(tsort(g), keys(m, :dep))
  V = vertices(g)
  lookup = Dict{String,Integer}()
  for v in V
    setindex!(lookup, v.index, v.label)
  end
  for i in 1:length(samplers)
    sampler = samplers[i]
    links = String[]
    for key in sampler.params
      append!(links, getlinks(V[lookup[key]], g, m))
    end
    sampler.links = intersect(m.links, links)
  end
  m
end


#################### MCMCModel Base Methods ####################

function Base.getindex(m::MCMCModel, key::String)
  m.nodes[key]
end

function Base.keys(m::MCMCModel, ntype::Symbol=:assigned, block::Integer=0)
  values = String[]
  if ntype == :all
    for key in keys(m.nodes)
      if isa(m[key], MCMCDepNode)
        values = [values, key, m[key].deps]
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
      if isa(m[key], MCMCDepNode)
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
      if isa(node, MCMCDepNode) && node.monitor
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

function labels{T<:String}(m::MCMCModel, nkeys::Vector{T})
  values = String[]
  for key in nkeys
    node = m[key]
    if isa(node, MCMCDepNode)
      if isa(node, VariateScalar)
        push!(values, key)
      elseif isa(node, VariateVector)
        for i in 1:length(node)
          push!(values, string(key, "[", i, "]"))
        end
      elseif isa(node, VariateMatrix)
        for j in 1:size(node, 2)
          for i in 1:size(node, 1)
            push!(values, string(key, "[", i, ",", j, "]"))
          end
        end
      else
        error("unsupported MCMCDepNode node")
      end
    else
      error("only MCMCDepNode nodes may be labeled")
    end
  end
  values
end

function setinits!{T<:String}(m::MCMCModel, inits::Dict{T,Any})
  for key in m.links
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
    isa(inputs[key], MCMCDepNode) && error("inputs must not be MCMCDepNode types")
    m.nodes[key] = deepcopy(inputs[key])
  end
  m.hasinputs = true
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
    append!(nkeys, m.samplers[b].links)
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

function unlist{T<:String}(m::MCMCModel, nkeys::Vector{T}, transform::Bool=false)
  f = transform ? link : identity
  N = map(key -> length(m[key]), nkeys)
  values = Array(VariateType, sum(N))
  i = 0
  for k in 1:length(nkeys)
    node = m[nkeys[k]]
    n = N[k]
    values[i+(1:n)] = f(node, node.data)
    i += n
  end
  values
end

function update!(m::MCMCModel, block::Integer=0)
  links = block > 0 ? m.samplers[block].links : m.links
  for key in links
    update!(m[key], m)
  end
  m
end
