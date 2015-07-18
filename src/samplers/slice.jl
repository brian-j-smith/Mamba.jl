#################### Slice Sampler ####################

#################### Types ####################

type SliceTune
  width::Vector{Float64}
end

type SliceVariate <: VectorVariate
  value::Vector{Float64}
  tune::SliceTune

  SliceVariate(x::Vector{Float64}, tune::SliceTune) = new(x, tune)
end

function SliceVariate(x::Vector{Float64}, tune=nothing)
  tune = SliceTune(
    Array(Float64, 0)
  )
  SliceVariate(x, tune)
end


#################### Sampler Constructor ####################

function Slice{T<:Real}(params::Vector{Symbol}, width::Vector{T},
           stype::Symbol=:multivar; transform::Bool=false)
  Sampler(params,
    quote
      tunepar = tune(model, block)
      x = unlist(model, block, tunepar["transform"])
      f = x -> logpdf!(model, x, block, tunepar["transform"])
      v = SliceVariate(x)
      slice!(v, tunepar["width"], f, tunepar["stype"])
      relist(model, v.value, block, tunepar["transform"])
    end,
    Dict("width" => Float64[width...], "stype" => stype,
         "transform" => transform)
  )
end


#################### Sampling Functions ####################

function slice!(v::SliceVariate, width::Vector{Float64}, logf::Function,
           stype::Symbol=:multivar)
  stype == :multivar ? slice_multi!(v, width, logf) :
  stype == :univar   ? slice_uni!(v, width, logf) :
    error("unsupported slice sampler type $stype")
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
