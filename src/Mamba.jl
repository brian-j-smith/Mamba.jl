using Distributions

module Mamba

  using Distributions


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
  import StatsBase: autocor, autocov, crosscov, describe, quantile, sem,
         StatsBase, summarystats


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
    monitor::Vector{Int}
    eval::Function
    sources::Vector{Symbol}
    targets::Vector{Symbol}
  end

  type Stochastic{T} <: Dependent{T}
    value::T
    symbol::Symbol
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
    tune::Dict
    targets::Vector{Symbol}
  end


  #################### Model Type ####################

  type Model
    nodes::Dict{Symbol,Any}
    dependents::Vector{Symbol}
    samplers::Vector{Sampler}
    state::Vector{Vector{VariateType}}
    iter::Integer
    burnin::Integer
    chain::Integer
    hasinputs::Bool
    hasinits::Bool
  end


  #################### Chains Type ####################

  immutable Chains
    value::Array{VariateType,3}
    names::Vector{String}
    range::Range{Int}
    model::Model
  end


  #################### Includes ####################

  include("distributions/constructors.jl")
  include("distributions/flat.jl")
  include("distributions/methods.jl")
  include("distributions/null.jl")

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
  include("output/mcse.jl")
  include("output/stats.jl")
  include("output/plot.jl")

  include("samplers/amm.jl")
  include("samplers/amwg.jl")
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

  export
    cm,
    inch,
    mm,
    pt,
    px


  #################### Deprecated ####################

  include("deprecated.jl")

end
