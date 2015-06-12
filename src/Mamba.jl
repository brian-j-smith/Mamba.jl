using Distributions

module Mamba

  using Distributions


  #################### User Add-on Packages and Functions ####################

  if isfile(joinpath(dirname(dirname(@__FILE__)),"usr","addons.jl"))
    include("../usr/addons.jl")
  end

  #################### Imports ####################

  import Base: Base, cor, dot
  import Base.LinAlg: Cholesky
  import Calculus: gradient
  import Compose: Context, context, cm, gridstack, inch, MeasureOrNumber, mm,
         pt, px
  import Distributions: Continuous, ContinuousUnivariateDistribution,
         Distribution, Distributions, gradlogpdf, insupport, logpdf, logpdf!,
         minimum, maximum, MatrixDistribution, MultivariateDistribution,
         PDiagMat, PDMat, quantile, rand, ScalMat, Truncated,
         UnivariateDistribution, ValueSupport
  import Gadfly: draw, Geom, Guide, Layer, layer, PDF, Plot, plot, PNG, PS,
         render, Scale, SVG, Theme
  import Graphs: AbstractGraph, add_edge!, add_vertex!, Edge, KeyVertex, graph,
         out_edges, out_neighbors, target, topological_sort_by_dfs, vertices
  import Showoff: showoff
  import StatsBase: autocor, autocov, counts, describe, predict, quantile, sem,
         StatsBase, summarystats

  include("distributions/pdmats2.jl")
  importall .PDMats2


  #################### Variate Types ####################

  typealias VariateType Float64

  abstract Variate{T<:Union(VariateType, Array{VariateType})}

  typealias UniVariate Variate{VariateType}
  typealias MultiVariate{N} Variate{Array{VariateType,N}}
  typealias VectorVariate MultiVariate{1}
  typealias MatrixVariate MultiVariate{2}


  #################### Distribution Types ####################

  typealias DistributionStruct Union(Distribution, Array{Distribution})


  #################### Dependent Types ####################

  abstract Dependent{T} <: Variate{T}

  type Logical{T} <: Dependent{T}
    value::T
    symbol::Symbol
    nlink::Integer
    monitor::Vector{Int}
    eval::Function
    sources::Vector{Symbol}
    targets::Vector{Symbol}
  end

  type Stochastic{T} <: Dependent{T}
    value::T
    symbol::Symbol
    nlink::Integer
    monitor::Vector{Int}
    eval::Function
    sources::Vector{Symbol}
    targets::Vector{Symbol}
    distr::DistributionStruct
  end


  #################### Sampler Type ####################

  type Sampler
    params::Vector{Symbol}
    eval::Function
    tune::Dict{String,Any}
    targets::Vector{Symbol}
  end


  #################### Model Type ####################

  type Model
    nodes::Dict{Symbol,Any}
    dependents::Vector{Symbol}
    samplers::Vector{Sampler}
    states::Vector{Vector{VariateType}}
    iter::Integer
    burnin::Integer
    chain::Integer
    hasinputs::Bool
    hasinits::Bool
  end


  #################### Chains Type ####################

  immutable Chains
    value::Array{VariateType,3}
    range::Range{Int}
    names::Vector{String}
    chains::Vector{Integer}
    model::Model
  end


  #################### Includes ####################

  include("distributions/flat.jl")
  include("distributions/mvnormal.jl")
  include("distributions/null.jl")
  include("distributions/constructors.jl")
  include("distributions/methods.jl")

  include("model/core.jl")
  include("model/dependent.jl")
  include("model/graph.jl")
  include("model/initialization.jl")
  include("model/mcmc.jl")
  include("model/simulation.jl")

  include("output/chains.jl")
  include("output/chainsummary.jl")
  include("output/gelmandiag.jl")
  include("output/gewekediag.jl")
  include("output/heideldiag.jl")
  include("output/mcse.jl")
  include("output/predict.jl")
  include("output/rafterydiag.jl")
  include("output/stats.jl")
  include("output/plot.jl")

  include("samplers/amm.jl")
  include("samplers/amwg.jl")
  include("samplers/dgs.jl")
  include("samplers/miss.jl")
  include("samplers/nuts.jl")
  include("samplers/sampler.jl")
  include("samplers/slice.jl")

  include("variate/core.jl")
  include("variate/numeric.jl")

  include("progress.jl")
  include("utils.jl")


  #################### Exports ####################

  export
    MatrixVariate,
    MultiVariate,
    UniVariate,
    Variate,
    VariateType,
    VectorVariate

  export
    Chains,
    Logical,
    Model,
    Sampler,
    Stochastic

  export
    BDiagNormal,
    Flat

  export
    @modelexpr,
    autocor,
    cor,
    describe,
    dic,
    draw,
    gelmandiag,
    gewekediag,
    gradlogpdf,
    gradlogpdf!,
    graph,
    graph2dot,
    heideldiag,
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
    predict,
    quantile,
    rafterydiag,
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
    DGS,
    MISS,
    nuts!,
    nutsepsilon,
    NUTS,
    NUTSVariate,
    slice!,
    Slice,
    SliceVariate

  export
    cm,
    inch,
    mm,
    pt,
    px


  #################### Deprecated ####################

  include("deprecated.jl")

end
