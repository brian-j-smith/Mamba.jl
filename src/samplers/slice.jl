#################### Slice Sampler ####################

#################### Types and Constructors ####################

type SliceTune <: SamplerTune
  width::Vector{Float64}

  function SliceTune(value::Vector{Float64}=Float64[])
    new(
      Array(Float64, 0)
    )
  end
end


typealias SliceVariate SamplerVariate{SliceTune}


#################### Sampler Constructor ####################

function Slice{T<:Real}(params::Vector{Symbol}, width::Vector{T},
                        stype::Symbol=:multivar; transform::Bool=false)
  width = Float64[width...]
  samplerfx = function(model::Model, block::Integer)
    v = SamplerVariate(model, block, transform)
    f = x -> logpdf!(model, x, block, transform)
    slice!(v, width, f, stype)
    relist(model, v, block, transform)
  end
  Sampler(params, samplerfx, SliceTune())
end


#################### Sampling Functions ####################

function slice!(v::SliceVariate, width::Vector{Float64}, logf::Function,
                stype::Symbol=:multivar)
  stype == :multivar ? slice_multi!(v, width, logf) :
  stype == :univar   ? slice_uni!(v, width, logf) :
    throw(ArgumentError("unsupported slice sampler type $stype"))
end

function slice_multi!(v::SliceVariate, width::Vector{Float64}, logf::Function)
  p0 = logf(v.value) + log(rand())

  n = length(v)
  lower = v - width .* rand(n)
  upper = lower + width

  x = width .* rand(n) + lower
  while logf(x) < p0
    for i in 1:n
      value = x[i]
      if value < v[i]
        lower[i] = value
      else
        upper[i] = value
      end
      x[i] = (upper[i] - lower[i]) * rand() + lower[i]
    end
  end
  v[:] = x
  v.tune.width = width

  v
end

function slice_uni!(v::SliceVariate, width::Vector{Float64}, logf::Function)
  logf0 = logf(v.value)
  for i in 1:length(v)
    p0 = logf0 + log(rand())

    lower = v[i] - width[i] * rand()
    upper = lower + width[i]

    x = v[i]
    v[i] = width[i] * rand() + lower
    while true
      logf0 = logf(v.value)
      logf0 < p0 || break
      value = v[i]
      if value < x
        lower = value
      else
        upper = value
      end
      v[i] = (upper - lower) * rand() + lower
    end
  end
  v.tune.width = width
  v
end
