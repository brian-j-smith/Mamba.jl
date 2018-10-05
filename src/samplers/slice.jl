#################### Slice Sampler ####################

#################### Types and Constructors ####################

const SliceForm = Union{Univariate, Multivariate}

struct SliceTune{F<:SliceForm} <: SamplerTune
  logf::Union{Function, Missing}
  width::Union{Float64, Vector{Float64}}

  SliceTune{F}() where {F<:SliceForm} = new{F}()

  SliceTune{F}(x::Vector, width) where {F<:SliceForm} =
    SliceTune{F}(x, width, missing)

  SliceTune{F}(x::Vector, width, logf::Function) where {F<:SliceForm} =
    SliceTune{F}(x, width, logf)

  SliceTune{F}(x::Vector, width::Real, logf::Union{Function, Missing}) where
    {F<:SliceForm} = new{F}(logf, Float64(width))

  SliceTune{F}(x::Vector, width::Vector, logf::Union{Function, Missing}) where
    {F<:SliceForm} = new{F}(logf, convert(Vector{Float64}, width))
end


const SliceUnivariate = SamplerVariate{SliceTune{Univariate}}
const SliceMultivariate = SamplerVariate{SliceTune{Multivariate}}

validate(v::SamplerVariate{SliceTune{F}}) where {F<:SliceForm} =
  validate(v, v.tune.width)

validate(v::SamplerVariate{SliceTune{F}}, width::Float64) where {F<:SliceForm} = v

function validate(v::SamplerVariate{SliceTune{F}}, width::Vector) where {F<:SliceForm}
  n = length(v)
  length(width) == n ||
    throw(ArgumentError("length(width) differs from variate length $n"))
  v
end


#################### Sampler Constructor ####################

function Slice(params::ElementOrVector{Symbol},
                width::ElementOrVector{T},
                ::Type{F}=Multivariate;
                transform::Bool=false) where {T<:Real, F<:SliceForm}
  samplerfx = function(model::Model, block::Integer)
    block = SamplingBlock(model, block, transform)
    v = SamplerVariate(block, width)
    sample!(v, x -> logpdf!(block, x))
    relist(block, v)
  end
  Sampler(params, samplerfx, SliceTune{F}())
end


#################### Sampling Functions ####################

sample!(v::Union{SliceUnivariate, SliceMultivariate}) = sample!(v, v.tune.logf)


function sample!(v::SliceUnivariate, logf::Function)
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


function sample!(v::SliceMultivariate, logf::Function)
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
