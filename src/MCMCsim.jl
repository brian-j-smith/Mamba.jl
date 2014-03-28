module MCMCsim

using Distributions


#################### Imports ####################

import Base: Base, cor
import Distributions: insupport, logpdf, quantile
import Graphs: AbstractGraph, add_edge!, add_vertex!, Edge, ExVertex, graph,
       out_edges, out_neighbors, target, topological_sort_by_dfs, vertices
import StatsBase: autocor, crosscov, describe, quantile, sem, StatsBase


#################### Variate Types ####################

typealias VariateType Float64

abstract Variate{T<:Union(VariateType, VecOrMat{VariateType})}

typealias VariateScalar Variate{VariateType}
typealias VariateVector Variate{Vector{VariateType}}
typealias VariateMatrix Variate{Matrix{VariateType}}

typealias Univariate VariateScalar
typealias Multivariate Union(VariateVector, VariateMatrix)


#################### Distribution Types ####################

typealias DistributionStruct Union(Distribution, VecOrMat{Distribution})


#################### MCMCParam Types ####################

abstract MCMCParam{T} <: Variate{T}

type MCMCLogical{T} <: MCMCParam{T}
  data::T
  monitor::Bool
  eval::Function
  deps::Vector{String}
end

type MCMCStochastic{T} <: MCMCParam{T}
  data::T
  monitor::Bool
  eval::Function
  deps::Vector{String}
  distr::DistributionStruct
  inits::Vector{T}
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
  hasdata::Bool
  hasinits::Bool
end


#################### MCMCChain Type ####################

type MCMCChain
  data::Array{VariateType,3}
  names::Vector{String}
  start::Integer
  thin::Integer
  model::MCMCModel
end


#################### Includes ####################

include("chain.jl")
include("densities.jl")
include("mcmc.jl")
include("model.jl")
include("nodes.jl")
include("sampler.jl")
include("utils.jl")
include("variate.jl")

include("samplers/amm.jl")
include("samplers/amwg.jl")
include("samplers/nuts.jl")
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
  blockkeys,
  blocktune,
  describe,
  dic,
  gelmandiag,
  gradient,
  gradient!,
  graph,
  graph2dot,
  hpd,
  initchain!,
  insupport,
  invlink,
  link,
  logpdf,
  logpdf!,
  mcmc,
  plot,
  relist,
  relist!,
  setdata!,
  setinits!,
  simulate!,
  unlist

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
