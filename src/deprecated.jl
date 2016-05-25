depwarn2(old::AbstractString, new::AbstractString, funcsym::Symbol) =
  Base.depwarn("$old is deprecated, use $new instead.", funcsym)

## Deprecated at 0.9.1

#################### Legacy Sampler Code ####################

export nutsepsilon

function nutsepsilon(v::NUTSVariate, logfgrad::Function)
  depwarn2("nutsepsilon(v::NUTSVariate, logfgrad::Function)",
           "NUTSVariate(x::AbstractVector{T<:Real}, logfgrad::Function)",
           :nutsepsilon)
  nutsepsilon(v.value, logfgrad)
end

function SamplerVariate(m::Model, block::Integer, transform::Bool=false)
  Base.depwarn("SamplerVariate(m::Model, block::Integer, transform::Bool) is deprecated.",
               :SamplerVariate)
  SamplerVariate(unlist(m, block, transform), m.samplers[block], m.iter)
end

function Slice{T<:Real}(params::ElementOrVector{Symbol},
                        width::ElementOrVector{T}, stype::Symbol;
                        transform::Bool=false)
  depwarn2("Slice(params, width, stype::Symbol)",
           "Slice(params, width, ::Type{F<:Union{Univariate, Multivariate}})",
           :Slice)
  F = stype == :univar   ? Univariate :
      stype == :multivar ? Multivariate :
        throw(ArgumentError("unsupported slice sampler type $stype"))
  Slice(params, width, F, transform=transform)
end
