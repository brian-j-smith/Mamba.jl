using Distributions

module Mamba

  using Distributions
  using Gadfly


  #################### Imports ####################

  import Base: Base, cor, dot
  import Base.LinAlg: Cholesky
  import Calculus: gradient
  import Compose: Context, context, cm, gridstack, inch, MeasureOrNumber, mm,
         pt, px
  import Distributions: Bernoulli, Categorical, Continuous,
         ContinuousUnivariateDistribution, Distribution, Distributions,
         gradlogpdf, insupport, logpdf, logpdf!, minimum, maximum,
         MatrixDistribution, MultivariateDistribution, PDiagMat, PDMat,
         quantile, rand, ScalMat, support, Truncated, UnivariateDistribution,
         ValueSupport
  import Gadfly: draw, Geom, Guide, Layer, layer, PDF, Plot, plot, PNG, PS,
         render, Scale, SVG, Theme
  import Graphs: AbstractGraph, add_edge!, add_vertex!, Edge, KeyVertex, graph,
         out_edges, out_neighbors, target, topological_sort_by_dfs, vertices
  import Showoff: showoff
  import StatsBase: autocor, autocov, countmap, counts, describe, predict,
         quantile, sem, StatsBase, summarystats

  include("distributions/pdmats2.jl")
  importall .PDMats2


  #################### Variate Types ####################

  abstract ScalarVariate <: Real
  abstract ArrayVariate{N} <: DenseArray{Float64,N}

  typealias AbstractVariate Union{ScalarVariate, ArrayVariate}
  typealias VectorVariate ArrayVariate{1}
  typealias MatrixVariate ArrayVariate{2}


  #################### Distribution Types ####################

  typealias DistributionStruct Union{Distribution,
                                     Array{UnivariateDistribution},
                                     Array{MultivariateDistribution}}


  #################### Dependent Types ####################

  type ScalarLogical <: ScalarVariate
    value::Float64
    symbol::Symbol
    linklength::Int
    monitor::Vector{Int}
    eval::Function
    sources::Vector{Symbol}
    targets::Vector{Symbol}
  end

  type ArrayLogical{N} <: ArrayVariate{N}
    value::Array{Float64,N}
    symbol::Symbol
    linklength::Int
    monitor::Vector{Int}
    eval::Function
    sources::Vector{Symbol}
    targets::Vector{Symbol}
  end

  type ScalarStochastic <: ScalarVariate
    value::Float64
    symbol::Symbol
    linklength::Int
    monitor::Vector{Int}
    eval::Function
    sources::Vector{Symbol}
    targets::Vector{Symbol}
    distr::UnivariateDistribution
  end

  type ArrayStochastic{N} <: ArrayVariate{N}
    value::Array{Float64,N}
    symbol::Symbol
    linklength::Int
    monitor::Vector{Int}
    eval::Function
    sources::Vector{Symbol}
    targets::Vector{Symbol}
    distr::DistributionStruct
  end

  typealias AbstractLogical Union{ScalarLogical, ArrayLogical}
  typealias AbstractStochastic Union{ScalarStochastic, ArrayStochastic}
  typealias AbstractDependent Union{AbstractLogical, AbstractStochastic}


  #################### Sampler Type ####################

  type Sampler
    params::Vector{Symbol}
    eval::Function
    tune::Dict{AbstractString,Any}
    targets::Vector{Symbol}
  end


  #################### Model Type ####################

  type Model
    nodes::Dict{Symbol,Any}
    dependents::Vector{Symbol}
    samplers::Vector{Sampler}
    states::Vector{Vector{Float64}}
    iter::Int
    burnin::Int
    chain::Int
    hasinputs::Bool
    hasinits::Bool
  end


  #################### Chains Type ####################

  abstract AbstractChains

  immutable Chains <: AbstractChains
    value::Array{Float64,3}
    range::Range{Int}
    names::Vector{AbstractString}
    chains::Vector{Int}
  end

  immutable ModelChains <: AbstractChains
    value::Array{Float64,3}
    range::Range{Int}
    names::Vector{AbstractString}
    chains::Vector{Int}
    model::Model
  end


  #################### Includes ####################

  include("progress.jl")
  include("utils.jl")
  include("variate.jl")

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
  include("output/modelchains.jl")
  include("output/modelstats.jl")
  include("output/rafterydiag.jl")
  include("output/stats.jl")
  include("output/plot.jl")

  include("samplers/amm.jl")
  include("samplers/amwg.jl")
  include("samplers/bmmg.jl")
  include("samplers/dgs.jl")
  include("samplers/miss.jl")
  include("samplers/nuts.jl")
  include("samplers/sampler.jl")
  include("samplers/slice.jl")


  #################### Exports ####################

  export
    AbstractChains,
    AbstractDependent,
    AbstractLogical,
    AbstractStochastic,
    AbstractVariate,
    ArrayLogical,
    ArrayStochastic,
    ArrayVariate,
    Chains,
    Logical,
    MatrixVariate,
    Model,
    ModelChains,
    Sampler,
    ScalarLogical,
    ScalarStochastic,
    ScalarVariate,
    Stochastic,
    VectorVariate

  export
    BDiagNormal,
    Flat

  export
    @modelexpr,
    autocor,
    changerate,
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
    bmmg!,
    BMMG,
    BMMGVariate,
    dgs!,
    DGS,
    DGSVariate,
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
