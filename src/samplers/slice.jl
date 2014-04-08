#################### Multivariate Slice Sampler ####################

#################### Types ####################

type TuneSlice
  width::Vector{Float64}
end

type VariateSlice <: VariateVector
  data::Vector{VariateType}
  tune::TuneSlice

  function VariateSlice{T<:Real}(x::Vector{T}, tune::TuneSlice)
    new(VariateType[x...], tune)
  end
end

function VariateSlice{T<:Real}(x::Vector{T}, tune=nothing)
  tune = TuneSlice(
    Array(Float64, 0)
  )
  VariateSlice(x, tune)
end


#################### Sampling Functions ####################

function slice{T<:Real}(x::Vector{T}, width::Vector{Float64}, logf::Function)
  slice!(VariateSlice(x), width, logf)
end

function slice!(v::VariateSlice, width::Vector{Float64}, logf::Function)
  p0 = logf(v.data) + log(rand())

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

#################### Sampling Functions ####################

function slicewg{T<:Real}(x::Vector{T}, width::Vector{Float64}, logf::Function)
  slicewg!(VariateSlice(x), width, logf)
end

function slicewg!(v::VariateSlice, width::Vector{Float64}, logf::Function)
  logf0 = logf(v.data)
  for i in 1:length(v)
    p0 = logf0 + log(rand())

    lower = v[i] - width[i] * rand()
    upper = lower + width[i]

    x = v[i]
    v[i] = width[i] * rand() + lower
    while true
      logf0 = logf(v.data)
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
