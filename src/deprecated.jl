## Deprecated at 0.8.1

function relist!{T<:Real}(m::Model, values::AbstractArray{T},
                 nodekeys::Vector{Symbol}, transform::Bool=false)
  msg = string("relist!{T<:Real}(m::Model, values::AbstractArray{T}, ",
               "nodekeys::Vector{Symbol}, transform::Bool=false) is deprecated")
  Base.depwarn(msg, :relist!)

  x = relist(m, values, nodekeys, transform)
  for key in nodekeys
    m[key].value = x[key]
  end
  update!(m)
end


export tune

@deprecate tune gettune


## Deprecated at 0.9.0

export simulate!

@deprecate simulate! sample!
