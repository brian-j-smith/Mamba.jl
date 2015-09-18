#################### Slice Simplex Sampler ####################

#################### Types ####################

type SliceSimplexTune
  scale::Float64
end

type SliceSimplexVariate <: VectorVariate
  value::Vector{Float64}
  tune::SliceSimplexTune

  SliceSimplexVariate(x::Vector{Float64}, tune::SliceSimplexTune) = new(x, tune)
end

function SliceSimplexVariate(x::Vector{Float64}, tune=nothing)
  tune = SliceSimplexTune(
    NaN
  )
  SliceSimplexVariate(x, tune)
end


#################### Sampler Constructor ####################

function SliceSimplex(params::Vector{Symbol}; scale::Real=1.0)
  Sampler(params,
    quote
      x = unlist(model, block)
      tunepar = tune(model, block)
      offset = 0
      for key in keys(model, :block, block)
        node = model[key]
        SliceSimplex_sub!(node.distr)
      end
      relist(model, x, block)
    end,
    Dict("scale" => scale)
  )
end

function SliceSimplex_sub!(D::Array{MultivariateDistribution})
  m = length(D)
  n = dims(D)[end]
  for i in 1:m
    inds = range(i, m, n) + offset
    v = SliceSimplexVariate(x[inds])
    f = function(y)
      x[inds] = y
      logpdf!(model, x, block)
    end
    slicesimplex!(v, f, scale=tunepar["scale"])
  end
  offset += m * n
end

function SliceSimplex_sub!(d::MultivariateDistribution)
  SliceSimplex_sub!(MultivariateDistribution[d])
end

function SliceSimplex_sub!(d)
  error("unsupported distribution type $(typeof(d)) in SliceSimplex")
end


#################### Sampling Functions ####################

function slicesimplex!(v::SliceSimplexVariate, logf::Function; scale::Real=1.0)
  d = Dirichlet(ones(v))

  insupport(d, v) || throw(ArgumentError("must supply a probability vector"))
  0 < scale <= 1 || throw(ArgumentError("scale must be > 0 and <= 1"))

  p0 = logf(v.value) + log(rand())
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
  vertices[:,2:end] += (1.0 - scale) * (vertices[:,1] .- vertices[:,2:end])
  vertices .+ x .- vertices * rand(Dirichlet(ones(x)))
end

function shrinksimplex(bx::AbstractVector{Float64}, bc::AbstractVector{Float64},
                       cx::AbstractVector{Float64}, cc::AbstractVector{Float64},
                       vertices::AbstractMatrix{Float64})
  for i in find(bc .< bx)
    inds = [1:(i-1); (i+1):size(vertices,2)]
    vertices[:,inds] += bc[i] * (vertices[:,i] .- vertices[:,inds])
    bc = vertices \ cc
  end
  vertices
end
