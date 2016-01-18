#################### Slice Simplex Sampler ####################

#################### Types and Constructors ####################

type SliceSimplexTune <: SamplerTune
  scale::Float64

  function SliceSimplexTune(value::Vector{Float64}=Float64[])
    new(
      NaN
    )
  end
end


typealias SliceSimplexVariate SamplerVariate{SliceSimplexTune}


#################### Sampler Constructor ####################

function SliceSimplex(params::ElementOrVector{Symbol}; scale::Real=1.0)
  params = asvec(params)
  samplerfx = function(model::Model, block::Integer)
    s = model.samplers[block]
    local node, x
    for key in params
      node = model[key]
      x = unlist(node)

      sim = function(inds::Range, logf::Function)
        v = SamplerVariate(x[inds], s, model.iter)
        slicesimplex!(v, logf, scale=scale)
      end

      logf = function(d::MultivariateDistribution, v::AbstractVector,
                      inds::Range)
        x[inds] = v
        relist!(model, x, key)
        logpdf(d, v) + logpdf(model, node.targets)
      end

      SliceSimplex_sub!(node.distr, sim, logf)
    end
    nothing
  end
  Sampler(params, samplerfx, SliceSimplexTune())
end


function SliceSimplex_sub!(d::MultivariateDistribution, sim::Function,
                           logf::Function)
  inds = 1:length(d)
  sim(inds, v -> logf(d, v, inds))
end

function SliceSimplex_sub!(D::Array{MultivariateDistribution}, sim::Function,
                           logf::Function)
  inds = 0:0
  for i in 1:length(D)
    d = D[i]
    inds = last(inds) + (1:length(d))
    sim(inds, v -> logf(d, v, inds))
  end
end

function SliceSimplex_sub!(d, sim::Function, logf::Function)
  throw(ArgumentError("unsupported distribution structure $(typeof(d))"))
end


#################### Sampling Functions ####################

function slicesimplex!(v::SliceSimplexVariate, logf::Function; scale::Real=1.0)
  0 < scale <= 1 || throw(ArgumentError("scale is not in (0, 1]"))

  p0 = logf(v.value) + log(rand())
  d = Dirichlet(ones(v))

  vertices = makefirstsimplex(v, scale)
  vb = vertices \ v
  xb = rand(d)
  x = vertices * xb
  while any(x .< 0.0) || any(x .> 1.0) || logf(x) < p0
    vertices = shrinksimplex(vb, xb, v, x, vertices)
    vb = vertices \ v
    xb = rand(d)
    x = vertices * xb
  end
  v[:] = x
  v.tune.scale = scale

  v
end


function makefirstsimplex(x::AbstractVector{Float64}, scale::Real)
  vertices = eye(length(x))
  vertices[:, 2:end] += (1.0 - scale) * (vertices[:, 1] .- vertices[:, 2:end])
  vertices .+ x .- vertices * rand(Dirichlet(ones(x)))
end


function shrinksimplex(bx::AbstractVector{Float64}, bc::AbstractVector{Float64},
                       cx::AbstractVector{Float64}, cc::AbstractVector{Float64},
                       vertices::AbstractMatrix{Float64})
  for i in find(bc .< bx)
    inds = [1:(i - 1); (i + 1):size(vertices, 2)]
    vertices[:, inds] += bc[i] * (vertices[:, i] .- vertices[:, inds])
    bc = vertices \ cc
  end
  vertices
end
