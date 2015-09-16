#################### Slice Simplex Sampler ####################

#################### Types ####################

type SliceSimplexTune
  width::Float64
end

type SliceSimplexVariate <: VectorVariate
  value::Vector{Float64}
  tune::SliceSimplexTune

  SliceSimplexVariate(x::Vector{Float64}, tune::SliceSimplexTune) = new(x, tune)
end

function SliceSimplexVariate(x::Vector{Float64}, tune=nothing)
  tune = SliceSimplexTune(1.0)
  SliceSimplexVariate(x, tune)
end


#################### Sampler Constructor ####################

function SliceSimplex{T<:Real}(params::Vector{Symbol}, 
                               width::T; transform::Bool=false)
  Sampler(params,
    quote
      tunepar = tune(model, block)
      x = unlist(model, block, tunepar["transform"])
      f = x -> logpdf!(model, x, block, tunepar["transform"])
      v = SliceSimplexVariate(x)
      slicesimplex!(v, tunepar["width"], f)
      relist(model, v.value, block, tunepar["transform"])
    end,
    ["width" => convert(Float64,width), "transform" => transform]
  )
end


#################### Sampling Functions ####################

function slicesimplex!(v::SliceSimplexVariate, width::Float64, logf::Function)
  p0 = logf(v.value) + log(rand())

  #simplex for generating theta
  vertices = makefirstsimplex(v.value,width)
  vb = vertices\v.value

  #candidate 
  xb = rand(Dirichlet(ones(length(v.value))))
  x = vertices * xb

  while logf(x) < p0 || !all(0.0 .< x .< 1)
    vertices = shrinksimplex(vb,xb,v.value,x,vertices)
    vb = vertices\v.value

    xb = rand(Dirichlet(ones(length(v.value))))
    x = vertices * xb
  end
  v[:] = x
  return v
end

#################### Helper ####################

function makefirstsimplex(kappa::Vector{Float64},width::Float64)
  vertices = diagm(ones(length(kappa)))

  if width < 1
    # get the coordinates of hte simplex shrunk twoards vertex 1
    vertices[:,2:end] = vertices[:,2:end] + 
      (1-width) * (vertices[:,1] .- vertices[:,2:end])

    # first, generate a point uniformly in appropriate size triangle
    rbaryc = rand(Dirichlet(ones(length(kappa))))

    # then translate vertices to make x correspond to random point
    vertices = vertices .+ kappa .- vertices * rbaryc
  end
  return vertices
end

function shrinksimplex(bx::Vector{Float64}, bc::Vector{Float64}, 
                       cx::Vector{Float64}, cc::Vector{Float64},
                       vertices::Matrix{Float64})
  for i in find(bc .< bx)
    vertices[:,[1:(i-1);(i+1):end]] = vertices[:, [1:(i-1);(i+1):end]] + 
      bc[i]*(vertices[:,i] .- vertices[:,[1:(i-1);(i+1):end]])

    bc = vertices\cc # this is solution to (vertices)*x = cc
  end
  return vertices
end
