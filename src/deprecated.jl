## Legacy type aliases

typealias VariateType Float64
typealias UniVariate ScalarVariate
typealias MultiVariate{N} ArrayVariate{N}


## Deprecated at 0.6.3

function link(d::AbstractDependent, x, transform::Bool=true)
  msg = string("link(d::$(typeof(d)), x, transform) is deprecated. ",
               "Use unlist(d, x, transform) instead.")
  Base.depwarn(msg, :link)
  unlist(d, x, transform)
end

function invlink(d::AbstractDependent, x, transform::Bool=true)
  msg = string("invlink(d::$(typeof(d)), x, transform) is deprecated. ",
               "Use relist(d, x, transform) instead.")
  Base.depwarn(msg, :invlink)
  relist(d, x, transform)
end
