#################### Multivariate Slice Sampler ####################

#################### Types ####################

type TuneSlice
  width::Vector{Float64}
end

type VariateSlice <: VariateVector
  data::Vector{VariateType}
  tune::TuneSlice

  VariateSlice(x::Vector, tune::TuneSlice) = new(VariateType[x...], tune)
end

function VariateSlice(x::Vector, tune=nothing)
  tune = TuneSlice(
    Array(Float64, 0)
  )
  VariateSlice(x, tune)
end


#################### Sampling Functions ####################

function slice(x::Vector, width::Vector{Float64}, logf::Function, args...)
  slice!(VariateSlice(x), width, logf, args...)
end

function slice!(v::VariateSlice, width::Vector{Float64}, logf::Function,
                args...)
  p0 = logf(v.data, args...) + log(rand())

  n = length(v)
  lower = v - width .* rand(n)
  upper = lower + width

  x = width .* rand(n) + lower
  while logf(x, args...) < p0
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

function slicewg(x::Vector, width::Vector{Float64}, logf::Function, args...)
  slicewg!(VariateSlice(x), width, logf, args...)
end

function slicewg!(v::VariateSlice, width::Vector{Float64}, logf::Function,
                  args...)
  logf0 = logf(v.data, args...)
  for i in 1:length(v)
    p0 = logf0 + log(rand())

    lower = v[i] - width[i] * rand()
    upper = lower + width[i]

    x = v[i]
    v[i] = width[i] * rand() + lower
    while true
      logf0 = logf(v.data, args...)
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
