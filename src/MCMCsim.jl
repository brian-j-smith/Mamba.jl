module MCMCsim

using Distributions


#################### Imports ####################

import Base: Base, cor
import Calculus: gradient
import Distributions: insupport, logpdf, PDiagMat, PDMat, quantile, ScalMat
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


#################### MCMCDepNode Types ####################

abstract MCMCDepNode{T} <: Variate{T}

type MCMCLogical{T} <: MCMCDepNode{T}
  data::T
  monitor::Bool
  eval::Function
  deps::Vector{String}
end

type MCMCStochastic{T} <: MCMCDepNode{T}
  data::T
  monitor::Bool
  eval::Function
  deps::Vector{String}
  distr::DistributionStruct
end


#################### MCMCSampler Type ####################

type MCMCSampler
  params::Vector{String}
  links::Vector{String}
  eval::Function
  tune::Dict
end


#################### MCMCModel Type ####################

type MCMCModel
  nodes::Dict{String,Any}
  links::Vector{String}
  samplers::Vector{MCMCSampler}
  iter::Integer
  burnin::Integer
  chain::Integer
  hasinputs::Bool
  hasinits::Bool
end


#################### MCMCChain Type ####################

type MCMCChains
  data::Array{VariateType,3}
  names::Vector{String}
  start::Integer
  thin::Integer
  model::MCMCModel
end


#################### Includes ####################

include("utils.jl")
include("variate.jl")

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


#################### Exports ####################

export
  Flat,
  MCMCChain,
  MCMCLogical,
  MCMCModel,
  MCMCSampler,
  MCMCStochastic

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
  link,
  logpdf,
  logpdf!,
  mcmc,
  plot,
  quantile,
  relist,
  relist!,
  setinits!,
  setinputs!,
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
  amm,
  amm!,
  amwg,
  amwg!,
  nuts,
  nuts!,
  nutseps,
  nutsfx!,
  slice,
  slice!,
  slicewg,
  slicewg!

end
