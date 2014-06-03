#################### Slice Sampler Types ####################

type SliceTune
  width::Vector{Float64}
end

type SliceVariate <: VectorVariate
  value::Vector{VariateType}
  tune::SliceTune

  SliceVariate(x::Vector{VariateType}, tune::SliceTune) = new(x, tune)
end

function SliceVariate(x::Vector{VariateType}, tune=nothing)
  tune = SliceTune(
    Array(Float64, 0)
  )
  SliceVariate(x, tune)
end


#################### Multivariate Slice Sampler ####################

#################### MCMCSampler Constructor ####################

function Slice{T<:Real}(params::Vector{Symbol}, width::Vector{T};
           transform::Bool=false)
  MCMCSampler(params,
    quote
      tunepar = tune(model, block)
      x = unlist(model, block, tunepar["transform"])
      f = x -> logpdf!(model, x, block, tunepar["transform"])
      v = SliceVariate(x)
      slice!(v, tunepar["width"], f)
      relist(model, v.value, block, tunepar["transform"])
    end,
    ["width" => Float64[width...], "transform" => transform]
  )
end


#################### Sampling Functions ####################

function slice!(v::SliceVariate, width::Vector{Float64}, logf::Function)
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


#################### Slice within Gibbs Sampler ####################

#################### Slice within Gibbs Sampler ####################

function SliceWG{T<:Real}(params::Vector{Symbol}, width::Vector{T};
           transform::Bool=false)
  MCMCSampler(params,
    quote
      tunepar = tune(model, block)
      x = unlist(model, block, tunepar["transform"])
      f = x -> logpdf!(model, x, block, tunepar["transform"])
      v = SliceVariate(x)
      slicewg!(v, tunepar["width"], f)
      relist(model, v.value, block, tunepar["transform"])
    end,
    ["width" => Float64[width...], "transform" => transform]
  )
end


#################### Sampling Functions ####################

function slicewg!(v::SliceVariate, width::Vector{Float64}, logf::Function)
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
