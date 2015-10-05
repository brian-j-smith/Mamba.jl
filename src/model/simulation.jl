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
  if block > 0
    sampler = m.samplers[block]
    params = sampler.params
    nkeys = union(sampler.params, sampler.targets)
  else
    params = keys(m, :block)
    nkeys = m.dependents
  end
  for key in nkeys
    value += logpdf(m[key], transform && in(key, params))
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
    isfinite(value) || return -Inf
  end
  for key in targets
    update!(m[key], m)
    value += logpdf(m[key], transform && in(key, params))
    isfinite(value) || return -Inf
  end
  value
end

function relist{T<:Real}(m::Model, values::AbstractArray{T}, block::Integer=0,
                         transform::Bool=false)
  relist(m, values, keys(m, :block, block), transform)
end

function relist{T<:Real}(m::Model, values::AbstractArray{T},
                         nkeys::Vector{Symbol}, transform::Bool=false)
  x = Dict{Symbol,Any}()
  N = length(values)
  offset = 0
  for key in nkeys
    value, n = relistlength(m[key], sub(values, (offset + 1):N), transform)
    x[key] = value
    offset += n
  end
  offset == length(values) ||
    throw(ErrorException("argument dimensions must match"))
  x
end

function relist!{T<:Real}(m::Model, values::AbstractArray{T}, block::Integer=0,
                          transform::Bool=false)
  relist!(m, convert(Array{Float64}, values), block, transform)
end

function relist!(m::Model, values::AbstractArray{Float64}, block::Integer=0,
                 transform::Bool=false)
  nkeys = keys(m, :block, block)
  x = relist(m, values, nkeys, transform)
  for key in nkeys
    m[key].value = x[key]
  end
  update!(m, block)
end

function relist!{T<:Real}(m::Model, values::AbstractArray{T},
                          nkeys::Vector{Symbol}, transform::Bool=false)
  relist!(m, convert(Array{Float64}, values), nkeys, transform)
end

function relist!(m::Model, values::AbstractArray{Float64},
                 nkeys::Vector{Symbol}, transform::Bool=false)
  x = relist(m, values, nkeys, transform)
  for key in nkeys
    m[key].value = x[key]
  end
  update!(m)
end

function simulate!(m::Model, block::Integer=0)
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

function unlist(m::Model, nkeys::Vector{Symbol}, transform::Bool=false)
  vcat(map(key -> unlist(m[key], transform), nkeys)...)
end

function update!(m::Model, block::Integer=0)
  nkeys = block > 0 ? m.samplers[block].targets : m.dependents
  for key in nkeys
    update!(m[key], m)
  end
  m
end
