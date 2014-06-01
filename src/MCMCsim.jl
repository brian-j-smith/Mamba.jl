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


#################### Variate Types ####################

typealias VariateType Float64

abstract AbstractVariate{T<:Union(VariateType, Array{VariateType})}

typealias AbstractVariateScalar AbstractVariate{VariateType}
typealias AbstractVariateVector AbstractVariate{Vector{VariateType}}
typealias AbstractVariateMatrix AbstractVariate{Matrix{VariateType}}
typealias AbstractVariateArray{N} AbstractVariate{Array{VariateType,N}}


#################### Distribution Types ####################

typealias DistributionStruct Union(Distribution, Array{Distribution})


#################### MCMCDependent Types ####################

abstract MCMCDependent{T} <: AbstractVariate{T}

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
  AbstractVariateArray,
  AbstractVariateMatrix,
  AbstractVariateScalar,
  AbstractVariateVector,
  MCMCChains,
  MCMCLogical,
  MCMCModel,
  MCMCSampler,
  MCMCStochastic,
  VariateType,
  Variate

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
  AMM,
  AMWG,
  MISS,
  NUTS,
  Slice,
  SliceWG,
  VariateAMM,
  VariateAMWG,
  VariateNUTS,
  VariateSlice

export
  amm!,
  amwg!,
  nuts!,
  nutsepsilon,
  slice!,
  slicewg!

end
