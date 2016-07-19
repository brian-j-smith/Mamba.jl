depwarn2(old::AbstractString, new::AbstractString, funcsym::Symbol) =
  Base.depwarn("$old is deprecated, use $new instead.", funcsym)
