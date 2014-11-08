function predict(c::Chains, key::Symbol)
  ismodelbased(c) || error("predict requires Chains from a Model fit")

  m = c.model
  node = m[key]
  isa(node, Stochastic) || error("predict is only defined for Stochastic nodes")

  nodenames = names(m, [key])

  sources = intersect(node.sources, keys(m, :stochastic))
  idx = indexin(names(m, sources), c.names)
  in(0, idx) && error("predict requires monitoring of nodes: ",
                      join(map(string, sources), ", "))

  iters, _, chains = size(c.value)
  value = Array(Float64, iters, length(nodenames), chains)
  for k in 1:chains
    for i in 1:iters
      relist!(m, vec(c.value[i,idx,k]), sources)
      v = isa(node.distr, Array) ? map(rand, node.distr) : rand(node.distr)
      value[i,:,k] = link(node, v, false)
    end
  end

  Chains(value, start=start(c.range), thin=step(c.range), names = nodenames)
end
