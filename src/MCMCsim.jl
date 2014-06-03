module MCMCsim

using Distributions




#################### Imports ####################

import Base: Base, cor, dot
import Base.LinAlg: Cholesky
import Calculus: gradient
import Distributions: Continuous, Distribution, insupport, logpdf, logpdf!,
       minimum, maximum, PDiagMat, PDMat, quantile, ScalMat, Truncated
import Graphs: AbstractGraph, add_edge!, add_vertex!, Edge, KeyVertex, graph,
       out_edges, out_neighbors, target, topological_sort_by_dfs, vertices
import StatsBase: autocor, autocov, crosscov, describe, quantile, sem,
       StatsBase, summarystats
import Gadfly: Plot, plot, Layer, layer, Geom, Guide, Theme, Scale
import Color: distinguishable_colors, LCHab



#################### Variate Types ####################

typealias VariateType Float64

abstract Variate{T<:Union(VariateType, Array{VariateType})}

typealias UniVariate Variate{VariateType}
typealias MultiVariate{N} Variate{Array{VariateType,N}}
typealias VectorVariate MultiVariate{1}
typealias MatrixVariate MultiVariate{2}


#################### Distribution Types ####################

typealias DistributionStruct Union(Distribution, Array{Distribution})


#################### MCMCDependent Types ####################

abstract MCMCDependent{T} <: Variate{T}

type MCMCLogical{T} <: MCMCDependent{T}
  value::T
  symbol::Symbol
  monitor::Vector{Bool}
  eval::Function
  sources::Vector{Symbol}
  targets::Vector{Symbol}
end

type MCMCStochastic{T} <: MCMCDependent{T}
  value::T
  symbol::Symbol
  monitor::Vector{Bool}
  eval::Function
  sources::Vector{Symbol}
  targets::Vector{Symbol}
  distr::DistributionStruct
end


#################### MCMCSampler Type ####################

type MCMCSampler
  params::Vector{Symbol}
  eval::Function
  tune::Dict
  targets::Vector{Symbol}
end


#################### MCMCModel Type ####################

type MCMCModel
  nodes::Dict{Symbol,Any}
  dependents::Vector{Symbol}
  samplers::Vector{MCMCSampler}
  iter::Integer
  burnin::Integer
  chain::Integer
  hasinputs::Bool
  hasinits::Bool
end


#################### MCMCChain Type ####################

immutable MCMCChains
  value::Array{VariateType,3}
  names::Vector{String}
  range::Range{Int}
  model::MCMCModel
end


#################### Includes ####################

include("distributions/constructors.jl")
include("distributions/flat.jl")
include("distributions/methods.jl")
include("distributions/null.jl")

include("model/graph.jl")
include("model/mcmc.jl")
include("model/model.jl")
include("model/nodes.jl")

include("output/chains.jl")
include("output/chainsummary.jl")
include("output/gelmandiag.jl")
include("output/gewekediag.jl")
include("output/mcse.jl")
include("output/stats.jl")
include("output/plot.jl")

include("samplers/amm.jl")
include("samplers/amwg.jl")
include("samplers/miss.jl")
include("samplers/nuts.jl")
include("samplers/sampler.jl")
include("samplers/slice.jl")

include("utils.jl")

include("variate/core.jl")
include("variate/numeric.jl")




#################### Exports ####################

export
  MatrixVariate,
  MultiVariate,
  UniVariate,
  Variate,
  VariateType,
  VectorVariate

export
  MCMCChains,
  MCMCLogical,
  MCMCModel,
  MCMCSampler,
  MCMCStochastic

export
  Flat

export
  @modelexpr,
  autocor,
  cor,
  describe,
  dic,
  gelmandiag,
  gewekediag,
  gradient,
  gradient!,
  graph,
  graph2dot,
  hpd,
  insupport,
  invlink,
  invlogit,
  link,
  logit,
  logpdf,
  logpdf!,
  mcmc,
  mcse,
  plot,
  quantile,
  relist,
  relist!,
  setinits!,
  setinputs!,
  setmonitor!,
  setsamplers!,
  simulate!,
  summarystats,
  tune,
  unlist,
  update!

export
  amm!,
  AMM,
  AMMVariate,
  amwg!,
  AMWG,
  AMWGVariate,
  MISS,
  nuts!,
  nutsepsilon,
  NUTS,
  NUTSVariate,
  slice!,
  Slice,
  slicewg!,
  SliceWG,
  SliceVariate

end
