#################### Slice Sampler ####################

#################### Types and Constructors ####################

type SliceTune <: SamplerTune
  width::Union{Real, Vector}

  function SliceTune(value::Vector{Float64}=Float64[])
    new(
      NaN
    )
  end
end


typealias SliceVariate SamplerVariate{SliceTune}


#################### Sampler Constructor ####################

function Slice{T<:Real}(params::Vector{Symbol}, width::ElementOrVector{T},
                        stype::Symbol=:multivar; transform::Bool=false)
  samplerfx = function(model::Model, block::Integer)
    v = SamplerVariate(model, block, transform)
    f = x -> logpdf!(model, x, block, transform)
    slice!(v, width, f, stype)
    relist(model, v, block, transform)
  end
  Sampler(params, samplerfx, SliceTune())
end


#################### Sampling Functions ####################

function slice!{T<:Real}(v::SliceVariate, width::ElementOrVector{T},
                         logf::Function, stype::Symbol=:multivar)
  v.tune.width = width
  stype == :multivar ? slice_multi!(v, logf) :
  stype == :univar   ? slice_uni!(v, logf) :
    throw(ArgumentError("unsupported slice sampler type $stype"))
end


function slice_multi!(v::SliceVariate, logf::Function)
  p0 = logf(v.value) + log(rand())

  n = length(v)
  lower = v - v.tune.width .* rand(n)
  upper = lower + v.tune.width

  x = v.tune.width .* rand(n) + lower
  while logf(x) < p0
    for i in 1:n
      value = x[i]
      if value < v[i]
        lower[i] = value
      else
        upper[i] = value
      end
      x[i] = rand(Uniform(lower[i], upper[i]))
    end
  end
  v[:] = x

  v
end


function slice_uni!(v::SliceVariate, logf::Function)
  logf0 = logf(v.value)

  n = length(v)
  lower = v - v.tune.width .* rand(n)
  upper = lower + v.tune.width

  for i in 1:n
    p0 = logf0 + log(rand())

    x = v[i]
    v[i] = rand(Uniform(lower[i], upper[i]))
    while true
      logf0 = logf(v.value)
      logf0 < p0 || break
      value = v[i]
      if value < x
        lower[i] = value
      else
        upper[i] = value
      end
      v[i] = rand(Uniform(lower[i], upper[i]))
    end
  end

  v
end
