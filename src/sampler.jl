#################### MCMCSampler Constructor ####################

function MCMCSampler(params::Vector, expr::Expr, tune::Dict=Dict())
  MCMCSampler(String[params...], String[], samplerfx(expr), tune)
end


#################### MCMCSampler Methods ####################

function Base.show(io::IO, s::MCMCSampler)
  print(io, "An object of type \"$(summary(s))\"\n")
  print(io, "Sampling Block Nodes:\n")
  show(io, s.params)
  print(io, "\n\n")
  show(io, s.eval.code)
  print(io, "\n")
end

function Base.showall(io::IO, s::MCMCSampler)
  show(io, s)
  print(io, "\nTuning Parameters:\n")
  show(io, s.tune)
  print(io, "\nDirectly Linked Nodes:\n")
  show(io, s.links)
  print(io, "\n")
end


#################### Adaptive Metropolis within Gibbs ####################

function SamplerAMM{T<:String}(params::Vector{T}, Sigma::Matrix;
                               adapt::Bool=false)
  MCMCSampler(params,
    quote
      keys = blockkeys(model, block)
      x = unlist(model, keys)
      tunepar = blocktune(model, block)
      v = VariateAMM(x, tunepar["sampler"])
      amm!(v, tunepar["Sigma"], logpdfxm!, model, block)
      tunepar["sampler"] = v.tune
      relist(model, v.data, keys)
    end,
    ["Sigma" => cholfact(Sigma), "adapt" => adapt, "sampler" => nothing]
  )
end


#################### Adaptive Metropolis within Gibbs ####################

function SamplerAMWG{T<:String}(params::Vector{T}, sigma::Vector;
                                adapt::Bool=false, batch::Integer=50,
                                target::Real=0.44)
  MCMCSampler(params,
    quote
      keys = blockkeys(model, block)
      x = unlist(model, keys)
      tunepar = blocktune(model, block)
      v = VariateAMWG(x, tunepar["sampler"])
      amwg!(v, tunepar["sigma"], logpdfxm!, model, block,
            adapt=tunepar["adapt"], batch=tunepar["batch"],
            target=tunepar["target"])
      tunepar["sampler"] = v.tune
      relist(model, v.data, keys)
    end,
    ["sigma" => sigma, "adapt" => adapt, "batch" => batch, "target" => target,
     "sampler" => nothing]
  )
end


#################### Multivariate Slice Sampler ####################

function SamplerSlice{T<:String}(params::Vector{T}, width::Vector)
  MCMCSampler(params,
    quote
      keys = blockkeys(model, block)
      x = unlist(model, keys)
      v = slice(x, blocktune(model, block)["width"], logpdfxm!, model, block)
      relist(model, v.data, keys)
    end,
    ["width" => width]
  )
end


#################### Slice within Gibbs Sampler ####################

function SamplerSliceWG{T<:String}(params::Vector{T}, width::Vector)
  MCMCSampler(params,
    quote
      keys = blockkeys(model, block)
      x = unlist(model, keys)
      v = slicewg(x, blocktune(model, block)["width"], logpdfxm!, model, block)
      relist(model, v.data, keys)
    end,
    ["width" => width]
  )
end


#################### Utility Functions ####################

function logpdfxm!(x::Vector, m::MCMCModel, block::Integer)
  keys = blockkeys(m, block)
  relist!(m, x, keys)
  if all(map(key -> insupport(m[key]), keys))
    update!(m, block)
    logpdf(m, block)
  else
    -Inf
  end
end

function samplerfx(expr::Expr)
  eval(Main, Expr(:function, :(model::MCMCModel, block::Integer), expr))
end
