#################### MCMCModel Constructor ####################

function MCMCModel(; iter::Integer=0, burnin::Integer=0, chain::Integer=1,
                   samplers::Vector=MCMCSampler[], params...)
  nodes = Dict{String,Any}()
  for (key, value) in params
    isa(value, MCMCNode) || error("params must be MCMCNode types")
    nodes[string(key)] = deepcopy(value)
  end
  m = MCMCModel(nodes, String[], samplers, iter, burnin, chain, false, false)
  g = graph(m)
  m.links = intersect(tsort(g), paramkeys(m))
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

function Base.keys(m::MCMCModel, monitoronly::Bool=false)
  if monitoronly
    result = String[]
    for key in keys(m.nodes)
      node = m[key]
      if isa(node, MCMCNode) && node.monitor
        push!(result, key)
      end
    end
  else
    result = collect(keys(m.nodes))
  end
  result
end

function Base.setindex!{T<:String}(m::MCMCModel, values::Dict, keys::Vector{T})
  for key in keys
    m[key][:] = values[key]
  end
end

function Base.setindex!{T<:String}(m::MCMCModel, value, keys::Vector{T})
  length(keys) == 1 || throw(BoundsError())
  m[keys[1]][:] = value
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

function blockkeys(m::MCMCModel, block::Integer=0)
  blocks = block > 0 ? block : 1:length(m.samplers)
  values = String[]
  for b in blocks
    append!(values, m.samplers[b].params)
  end
  unique(values)
end

function blocktune(m::MCMCModel, block::Integer=0)
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

function inputkeys(m::MCMCModel)
  setdiff(nodekeys(m), paramkeys(m))
end

function labels{T<:String}(m::MCMCModel, keys::Vector{T})
  values = String[]
  for key in keys
    node = m[key]
    if isa(node, MCMCNode)
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
        error("unsupported MCMCNode node")
      end
    else
      error("only MCMCNode nodes may be labeled")
    end
  end
  values
end

function nodekeys(m::MCMCModel)
  result = String[]
  for key in paramkeys(m)
    result = [result, key, m[key].deps]
  end
  unique(result)
end

function paramkeys(m::MCMCModel)
  result = String[]
  for key in keys(m)
    if isa(m[key], MCMCNode)
      push!(result, key)
    end
  end
  result
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
  for key in inputkeys(m)
    isa(inputs[key], MCMCNode) && error("inputs must not be MCMCNode types")
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

#################### MCMCModel Simulation Methods ####################

function gradient(m::MCMCModel, block::Integer=0, transform::Bool=false,
                  dtype::Symbol=:central)
  x0 = unlist(m, block, transform)
  value = gradient!(m, x0, block, transform, dtype)
  relist!(m, x0, block, transform)
  value
end

function gradient(m::MCMCModel, x::Vector, block::Integer=0,
                  transform::Bool=false, dtype::Symbol=:central)
  x0 = unlist(m, block)
  value = gradient!(m, x, block, transform, dtype)
  relist!(m, x0, block)
  value
end

function gradient!(m::MCMCModel, x::Vector, block::Integer=0,
                   transform::Bool=false, dtype::Symbol=:central)
  f = x -> logpdf!(m, x, block, transform)
  gradient(f, x, dtype)
end

function logpdf(m::MCMCModel, block::Integer=0, transform::Bool=false)
  blocks = block > 0 ? block : 1:length(m.samplers)
  keys = String[]
  for b in blocks
    append!(keys, m.samplers[b].params)
    append!(keys, m.samplers[b].links)
  end
  mapreduce(key -> logpdf(m[key], transform), +, unique(keys))
end

function logpdf(m::MCMCModel, x::Vector, block::Integer=0,
                transform::Bool=false)
  x0 = unlist(m, block)
  value = logpdf!(m, x, block, transform)
  relist!(m, x0, block)
  value
end

function logpdf!(m::MCMCModel, x::Vector, block::Integer=0,
                 transform::Bool=false)
  keys = blockkeys(m, block)
  m[keys] = relist(m, x, keys, transform)
  if all(map(key -> insupport(m[key]), keys))
    update!(m, block)
    logpdf(m, block, transform)
  else
    -Inf
  end
end

function relist(m::MCMCModel, values::Vector, block::Integer=0,
                transform::Bool=false)
  relist(m, values, blockkeys(m, block), transform)
end

function relist{T<:String}(m::MCMCModel, values::Vector, keys::Vector{T},
                           transform::Bool=false)
  f =  transform ? invlink : identity
  x = Dict{String,Any}()
  j = 0
  for key in keys
    node = m[key]
    n = length(node)
    x[key] = f(node, values[j+(1:n)])
    j += n
  end
  j == length(values) || throw(ErrorException("argument dimensions must match"))
  x
end

function relist!(m::MCMCModel, values::Vector, block::Integer=0,
                 transform::Bool=false)
  keys = blockkeys(m, block)
  m[keys] = relist(m, values, keys, transform)
  update!(m, block)
end

function relist!{T<:String}(m::MCMCModel, values::Vector, keys::Vector{T},
                            transform::Bool=false)
  m[keys] = relist(m, values, keys, transform)
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
  unlist(m, blockkeys(m, block), transform)
end

function unlist{T<:String}(m::MCMCModel, keys::Vector{T}, transform::Bool=false)
  f = transform ? link : identity
  N = map(key -> length(m[key]), keys)
  values = Array(VariateType, sum(N))
  i = 0
  for k in 1:length(keys)
    node = m[keys[k]]
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
