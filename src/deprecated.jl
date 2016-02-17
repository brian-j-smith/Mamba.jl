## Deprecated at 0.8.1

function relist!{T<:Real}(m::Model, values::AbstractArray{T},
                 nodekeys::Vector{Symbol}, transform::Bool=false)
  msg = string("relist!{T<:Real}(m::Model, values::AbstractArray{T}, ",
               "nodekeys::Vector{Symbol}, transform::Bool=false) is deprecated")
  Base.depwarn(msg, :relist!)

  x = relist(m, values, nodekeys, transform)
  for key in nodekeys
    m[key].value = x[key]
  end
  update!(m)
end


export tune

@deprecate tune gettune


## Deprecated at 0.9.0

export simulate!

@deprecate simulate! sample!


## Legacy Sampler Code

export amm!, amwg!, bhmc!, bmc3!, bmg!, dgs!, hmc!, mala!, nuts!, rwm!, slice!,
       SliceVariate, slicesimplex!

depwarn2(old::AbstractString, new::AbstractString, funcsym::Symbol) =
  Base.depwarn("$old is deprecated, use $new instead.", funcsym)

depwarnsample(old::AbstractString) = depwarn2(old, "sample!(v)", symbol(old))

## AMM

function AMMTune(x::Vector)
  depwarn2("AMMVariate(x)", "AMMVariate(x, Sigma, logf)", :AMMVariate)
  n = length(x)
  AMMTune(x, fill(NaN, n, n))
end

function amm!(v::AMMVariate, SigmaF::Cholesky{Float64}, logf::Function;
              adapt::Bool=true)
  depwarnsample("amm!")
  if v.tune.m == 0
    v.tune.SigmaL = SigmaF[:L]
  end
  sample!(v, logf, adapt=adapt)
end

#AMWG

function AMWGTune(x::Vector)
  depwarn2("AMWGVariate(x)", "AMWGVariate(x, sigma, logf)", :AMWGVariate)
  AMWGTune(x, NaN)
end

function amwg!{T<:Real}(v::AMWGVariate, sigma::Vector{T}, logf::Function;
                        adapt::Bool=true, batchsize::Integer=50,
                        target::Real=0.44)
  depwarnsample("amwg!")
  if v.tune.m == 0
    v.tune.sigma = adapt ? copy(sigma) : sigma
  end
  v.tune.batchsize = batchsize
  v.tune.target = target
  sample!(v, logf, adapt=adapt)
end

## BHMC

function BHMCTune(x)
  depwarn2("BHMCVariate(x)", "BHMCVariate(x, traveltime, logf)", :BHMCVariate)
  BHMCTune(x, NaN)
end

function bhmc!(v::BHMCVariate, traveltime::Real, logf::Function)
  depwarnsample("bhmc!")
  v.tune.traveltime = traveltime
  sample!(v, logf)
end

## BHMC3

function BMC3(params::ElementOrVector{Symbol}, indexset::Vector{Vector{Int}})
  depwarn2("BMC3(params, indexset)", "BMC3(params)", :BMC3)
  samplerfx = function(model::Model, block::Integer)
    block = SamplingBlock(model, block)
    v = SamplerVariate(block, 0)
    bmc3!(v, indexset, x -> logpdf!(block, x))
    relist(block, v)
  end
  Sampler(params, samplerfx, BMC3Tune())
end

function bmc3!(v::BMC3Variate, logf::Function; k::Integer=1)
  depwarnsample("bmc3!")
  v.tune.k = k
  sample!(v, logf)
end

function bmc3!(v::BMC3Variate, indexset::Vector{Vector{Int}}, logf::Function)
  depwarnsample("bmc3!")
  v.tune.indexset = indexset
  x = v[:]
  idx = rand(v.tune.indexset)
  x[idx] = 1.0 - v[idx]
  if rand() < exp(logf(x) - logf(v.value))
    v[:] = x
  end
  v
end

## BMG

function bmg!(v::BMGVariate, logf::Function; k::Integer=1)
  depwarnsample("bmg!")
  k <= length(v) || throw(ArgumentError("k is greater than length(v)"))
  v.tune.k = k
  sample!(v, logf)
end

## DGS

function dgs!{T<:Real}(v::DGSVariate, support::Matrix{T}, logf::Function)
  depwarnsample("dgs!")
  v.tune.support = support
  sample!(v, x -> exp(logf(x)))
end

function dgs!{T<:Real}(v::DGSVariate, support::Matrix{T},
                       probs::Vector{Float64})
  depwarnsample("dgs!")
  v.tune.support = support
  sample!(v, probs)
end

## HMC

function HMCTune(x)
  depwarn2("HMCVariate(x)",
           "HMCVariate(x, epsilon, L, logfgrad) or" *
           " HMCVariate(x, epsilon, L, Sigma, logfgrad)",
           :HMCVariate)
  HMCTune(x, NaN, 0)
end

function hmc!(v::HMCVariate, epsilon::Real, L::Integer, logfgrad::Function)
  depwarnsample("hmc!")
  v.tune.epsilon = epsilon
  v.tune.L = L
  v.tune.SigmaL = I
  sample!(v, logfgrad)
end

function hmc!(v::HMCVariate, epsilon::Real, L::Integer,
              SigmaF::Cholesky{Float64}, logfgrad::Function)
  depwarnsample("hmc!")
  v.tune.epsilon = epsilon
  v.tune.L = L
  v.tune.SigmaL = SigmaF[:L]
  sample!(v, logfgrad)
end

## MALA

function MALATune(x)
  depwarn2("MALAVariate(x)",
           "MALAVariate(x, scale, logfgrad) or" *
           " MALAVariate(x, scale, Sigma, logfgrad)",
           :MALAVariate)
  MALATune(x, NaN)
end

function mala!(v::MALAVariate, scale::Real, logfgrad::Function)
  depwarnsample("mala!")
  v.tune.scale = scale
  v.tune.SigmaL = I
  sample!(v, logfgrad)
end

function mala!(v::MALAVariate, scale::Real, SigmaF::Cholesky{Float64},
               logfgrad::Function)
  depwarnsample("mala!")
  v.tune.scale = scale
  v.tune.SigmaL = SigmaF[:L]
  sample!(v, logfgrad)
end

## NUTS

function NUTSTune(x)
  depwarn2("NUTSVariate(x)",
           "NUTSVariate(x, epsilon, logfgrad) or NUTSVariate(x, logfgrad)",
           :NUTSVariate)
  NUTSTune(x, NaN)
end

function nuts!(v::NUTSVariate, epsilon::Real, logfgrad::Function;
               adapt::Bool=false, target::Real=0.6)
  depwarnsample("nuts!")
  if v.tune.m == 0
    v.tune.epsilon = epsilon
  end
  v.tune.target = target
  sample!(v, logfgrad, adapt=adapt)
end

## RWM

function RWMTune(x)
  depwarn2("RWMVariate(x)", "RWMVariate(x, scale, logf)", :RWMVariate)
  RWMTune(x, NaN)
end

function rwm!{T<:Real}(v::RWMVariate, scale::ElementOrVector{T},
                       logf::Function; proposal::SymDistributionType=Normal)
  depwarnsample("rwm!")
  v.tune.scale = scale
  v.tune.proposal = proposal
  sample!(v, logf)
end

## Slice

function SliceVariate{T<:Real}(x::AbstractVector{T})
  depwarn2("SliceVariate(x)",
           "SliceUnivariate(x, width, logf) or" *
           " SliceMultivariate(x, width, logf)",
           :SliceVariate)
  SamplerVariate{SliceTune{SliceForm}}(x, NaN)
end

function slice!{T<:Real}(v::SamplerVariate{SliceTune{SliceForm}},
                         width::ElementOrVector{T}, logf::Function,
                         stype::Symbol=:multivar)
  depwarnsample("slice!")
  v.tune.width = width
  S = stype == :univar   ? Univariate :
      stype == :multivar ? Multivariate :
        throw(ArgumentError("unsupported slice sampler type $stype"))
  sample!(SamplerVariate{SliceTune{S}}(v, width), logf)
  v
end

## SliceSimplex

function slicesimplex!(v::SliceSimplexVariate, logf::Function; scale::Real=1.0)
  depwarnsample("slicesimplex!")
  0 < scale <= 1 || throw(ArgumentError("scale is not in (0, 1]"))
  v.tune.scale = scale
  sample!(v, logf)
end
