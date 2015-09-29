#################### ModelChains ####################

#################### Constructors ####################

function ModelChains(c::Chains, m::Model)
  ModelChains(c.value, c.range, c.names, c.chains, m)
end


#################### Conversions ####################

Base.convert(::Type{Chains}, mc::ModelChains) =
  Chains(mc.value, mc.range, mc.names, mc.chains)


#################### Indexing ####################

function Base.getindex(mc::ModelChains, window, names, chains)
  c = getindex(convert(Chains, mc), window, names, chains)
  ModelChains(c, mc.model)
end


#################### Auxilliary Functions ####################

function link(c::ModelChains, transform::Bool=true)
  if transform
    cc = copy(c.value)
    inds_queue = 1:length(c.names)
    for key in intersect(keys(c.model, :monitor), keys(c.model, :stochastic))
      node = c.model[key]
      inds = findin(c.names, names(node))
      if length(inds) > 0
        f(x) = unlist(node, link(node, relist(node, x)))
        cc[:,inds,:] = mapslices(f, cc[:,inds,:], 2)
        inds_queue = setdiff(inds_queue, inds)
      end
    end
    for j in inds_queue
      x = cc[:,j,:]
      if minimum(x) > 0.0
        cc[:,j,:] = maximum(x) < 1.0 ? logit(x) : log(x)
      end
    end
  else
    cc = c.value
  end
  cc
end
