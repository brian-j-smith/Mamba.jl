module MCMCsim

using Distributions


#################### Imports ####################

import Base: Base, cor, dot
import Calculus: gradient
import Distributions: Continuous, Distribution, insupport, logpdf, minimum,
       maximum, PDiagMat, PDMat, quantile, ScalMat, Truncated
import Graphs: AbstractGraph, add_edge!, add_vertex!, Edge, ExVertex, graph,
       out_edges, out_neighbors, target, topological_sort_by_dfs, vertices
import StatsBase: autocor, crosscov, describe, quantile, sem, StatsBase,
       summarystats


#################### Variate Types ####################

typealias VariateType Float64

abstract Variate{T<:Union(VariateType, VecOrMat{VariateType})}

typealias VariateScalar Variate{VariateType}
typealias VariateVector Variate{Vector{VariateType}}
typealias VariateMatrix Variate{Matrix{VariateType}}
typealias VariateVecOrMat Union(VariateVector, VariateMatrix)


#################### Distribution Types ####################

typealias DistributionStruct Union(Distribution, VecOrMat{Distribution})


#################### MCMCDependent Types ####################

abstract MCMCDependent{T} <: Variate{T}

type MCMCLogical{T} <: MCMCDependent{T}
  value::T
  name::String
  monitor::Vector{Bool}
  eval::Function
  sources::Vector{String}
  targets::Vector{String}
end

type MCMCStochastic{T} <: MCMCDependent{T}
  value::T
  name::String
  monitor::Vector{Bool}
  eval::Function
  sources::Vector{String}
  targets::Vector{String}
  distr::DistributionStruct
end


#################### MCMCSampler Type ####################

type MCMCSampler
  params::Vector{String}
  eval::Function
  tune::Dict
  targets::Vector{String}
end


#################### MCMCModel Type ####################

type MCMCModel
  nodes::Dict{String,Any}
  samplers::Vector{MCMCSampler}
  targets::Vector{String}
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
  start::Integer
  thin::Integer
  model::MCMCModel
end


#################### Includes ####################

include("distributions/constructors.jl")
include("distributions/methods.jl")

include("model/graph.jl")
include("model/mcmc.jl")
include("model/model.jl")
include("model/nodes.jl")

include("output/chains.jl")
include("output/chainsummary.jl")
include("output/gelmandiag.jl")
include("output/stats.jl")

include("samplers/amm.jl")
include("samplers/amwg.jl")
include("samplers/nuts.jl")
include("samplers/sampler.jl")
include("samplers/slice.jl")

include("utils.jl")

include("variate/core.jl")
include("variate/numeric.jl")


#################### Exports ####################

export
  Flat

export
  MCMCChains,
  MCMCLogical,
  MCMCModel,
  MCMCSampler,
  MCMCStochastic,
  VariateType,
  VariateScalar,
  VariateVector,
  VariateMatrix,
  VariateVecOrMat

export
  autocor,
  cor,
  describe,
  dic,
  gelmandiag,
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
  SamplerAMM,
  SamplerAMWG,
  SamplerNUTS,
  SamplerSlice,
  SamplerSliceWG,
  VariateAMM,
  VariateAMWG,
  VariateNUTS,
  VariateSlice

export
  amm!,
  amwg!,
  nuts!,
  nutseps,
  nutsfx!,
  slice!,
  slicewg!

end
