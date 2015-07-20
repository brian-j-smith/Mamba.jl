#################### ModelChains Constructor ####################

function ModelChains(c::Chains, m::Model)
  ModelChains(c.value, c.range, c.names, c.chains, m)
end


#################### ModelChains Indexing ####################

function Base.getindex(mc::ModelChains, window, names, chains)
  inds1 = window2inds(mc, window)
  inds2 = names2inds(mc, names)
  c = Chains(mc.value[inds1, inds2, chains],
             start = first(mc.range) + (first(inds1) - 1) * step(mc.range),
             thin = step(inds1) * step(mc.range), names = mc.names[inds2],
             chains = mc.chains[chains])
  ModelChains(c, mc.model)
end
