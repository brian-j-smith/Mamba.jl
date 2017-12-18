depwarn2(old::AbstractString, new::AbstractString, funcsym::Symbol) =
  Base.depwarn("$old is deprecated, use $new instead.", funcsym)

## Deprecated at 0.11.2

@deprecate logit(x::AbstractArray) logit.(x)
@deprecate invlogit(x::AbstractArray) invlogit.(x)
