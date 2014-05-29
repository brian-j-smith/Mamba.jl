module MCMCsim

using Distributions
using Gadfly
using KernelDensity


#################### Imports ####################

import Base: Base, cor, dot
import Base.LinAlg: Cholesky
import Calculus: gradient
import Distributions: Continuous, Distribution, insupport, logpdf, logpdf!,
       minimum, maximum, PDiagMat, PDMat, quantile, ScalMat, Truncated
import Graphs: AbstractGraph, add_edge!, add_vertex!, Edge, ExVertex, graph,
       out_edges, out_neighbors, target, topological_sort_by_dfs, vertices
import StatsBase: autocor, autocov, crosscov, describe, quantile, sem,
       StatsBase, summarystats
import Gadfly: Plot, plot
import Color: distinguishable_colors
import KernelDensity: kde


#################### Variate Types ####################

typealias VariateType Float64

abstract Variate{T<:Union(VariateType, Array{VariateType})}

typealias VariateScalar Variate{VariateType}
typealias VariateVector Variate{Vector{VariateType}}
typealias VariateMatrix Variate{Matrix{VariateType}}
typealias VariateArray{N} Variate{Array{VariateType,N}}


#################### Distribution Types ####################

typealias DistributionStruct Union(Distribution, Array{Distribution})


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
  sources::Vector{String}
  targets::Vector{String}
end


#################### MCMCModel Type ####################

type MCMCModel
  nodes::Dict{String,Any}
  dependents::Vector{String}
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

include("plots/chain_plot.jl")


#################### Exports ####################

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
  VariateArray

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
  
export
  chain_plot

end
