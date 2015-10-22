#################### Model Simulation ####################

function gradlogpdf(m::Model, block::Integer=0, transform::Bool=false;
                    dtype::Symbol=:forward)
  x0 = unlist(m, block, transform)
  value = gradlogpdf!(m, x0, block, transform, dtype=dtype)
  relist!(m, x0, block, transform)
  value
end

function gradlogpdf{T<:Real}(m::Model, x::AbstractArray{T}, block::Integer=0,
                             transform::Bool=false; dtype::Symbol=:forward)
  x0 = unlist(m, block)
  value = gradlogpdf!(m, x, block, transform, dtype=dtype)
  relist!(m, x0, block)
  value
end

function gradlogpdf!{T<:Real}(m::Model, x::AbstractArray{T}, block::Integer=0,
                              transform::Bool=false; dtype::Symbol=:forward)
  f = x -> logpdf!(m, x, block, transform)
  gradient(f, x, dtype)
end

function logpdf(m::Model, block::Integer=0, transform::Bool=false)
  value = 0.0
  if block != 0
    sampler = m.samplers[block]
    params = sampler.params
    nodekeys = union(sampler.params, sampler.targets)
  else
    params = keys(m, :block)
    nodekeys = keys(m, :stochastic)
  end
  for key in nodekeys
    value += logpdf(m[key], transform && key in params)
    isfinite(value) || return -Inf
  end
  value
end

function logpdf{T<:Real}(m::Model, x::AbstractArray{T}, block::Integer=0,
                         transform::Bool=false)
  x0 = unlist(m, block)
  value = logpdf!(m, x, block, transform)
  relist!(m, x0, block)
  value
end

function logpdf!{T<:Real}(m::Model, x::AbstractArray{T}, block::Integer=0,
                          transform::Bool=false)
  value = 0.0
  if block != 0
    sampler = m.samplers[block]
    params = sampler.params
    targets = sampler.targets
  else
    params = keys(m, :block)
    targets = keys(m, :stochastic)
  end
  m[params] = relist(m, x, params, transform)
  for key in setdiff(params, targets)
    value += logpdf(m[key], transform)
    isfinite(value) || return -Inf
  end
  for key in targets
    update!(m[key], m)
    value += logpdf(m[key], transform && key in params)
    isfinite(value) || return -Inf
  end
  value
end

function relist{T<:Real}(m::Model, values::AbstractArray{T}, block::Integer=0,
                         transform::Bool=false)
  relist(m, values, keys(m, :block, block), transform)
end

function relist{T<:Real}(m::Model, values::AbstractArray{T},
                         nodekeys::Vector{Symbol}, transform::Bool=false)
  x = Dict{Symbol,Any}()
  N = length(values)
  offset = 0
  for key in nodekeys
    value, n = relistlength(m[key], sub(values, (offset + 1):N), transform)
    x[key] = value
    offset += n
  end
  offset == length(values) ||
    throw(ArgumentError("incompatible number of values to put in nodes"))
  x
end

function relist!{T<:Real}(m::Model, values::AbstractArray{T}, block::Integer=0,
                          transform::Bool=false)
  relist!(m, convert(Array{Float64}, values), block, transform)
end

function relist!(m::Model, values::AbstractArray{Float64}, block::Integer=0,
                 transform::Bool=false)
  nodekeys = keys(m, :block, block)
  x = relist(m, values, nodekeys, transform)
  for key in nodekeys
    m[key].value = x[key]
  end
  update!(m, block)
end

function relist!{T<:Real}(m::Model, values::AbstractArray{T},
                          nodekeys::Vector{Symbol}, transform::Bool=false)
  relist!(m, convert(Array{Float64}, values), nodekeys, transform)
end

function relist!(m::Model, values::AbstractArray{Float64},
                 nodekeys::Vector{Symbol}, transform::Bool=false)
  x = relist(m, values, nodekeys, transform)
  for key in nodekeys
    m[key].value = x[key]
  end
  update!(m)
end

function simulate!(m::Model, block::Integer=0)
  if block != 0
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

function tune(m::Model, block::Integer=0)
  if block != 0
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

function unlist(m::Model, block::Integer=0, transform::Bool=false)
  unlist(m, keys(m, :block, block), transform)
end

function unlist(m::Model, monitoronly::Bool)
  f = function(key)
    node = m[key]
    lvalue = unlist(node)
    monitoronly ? lvalue[node.monitor] : lvalue
  end
  vcat(map(f, keys(m, :dependent))...)
end

function unlist(m::Model, nodekeys::Vector{Symbol}, transform::Bool=false)
  vcat(map(key -> unlist(m[key], transform), nodekeys)...)
end

function update!(m::Model, block::Integer=0)
  nodekeys = block != 0 ? m.samplers[block].targets : keys(m, :dependent)
  update!(m, nodekeys)
end

function update!(m::Model, nodekeys::Vector{Symbol})
  for key in nodekeys
    update!(m[key], m)
  end
  m
end
