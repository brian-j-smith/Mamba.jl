function slicewg!(v::SliceVariate, width::Vector{Float64}, logf::Function)
  warn("slicewg! is deprecated, use slice! with :univar option instead.")
  slice!(v, width, logf, :univar)
end

function SliceWG{T<:Real}(params::Vector{Symbol}, width::Vector{T};
           transform::Bool=false)
  warn("SliceWG is deprecated, use Slice with :univar option instead.")
  Slice(params, width, :univar, transform=transform)
end

export
  slicewg!,
  SliceWG


@deprecate MCMCChains Chains
@deprecate MCMCLogical Logical
@deprecate MCMCModel Model
@deprecate MCMCSampler Sampler
@deprecate MCMCStochastic Stochastic

export
  MCMCChains,
  MCMCLogical,
  MCMCModel,
  MCMCSampler,
  MCMCStochastic
