function predict(c::Chains, key::Symbol)
  ismodelbased(c) || error("predict requires Chains from a Model fit")

  m = c.model
  node = m[key]
  nodenames = names(m, [key])

  idx = indexin(names(m, node.sources), c.names)
  in(0, idx) && error("predict requires monitoring of nodes: ",
                      join(map(string, node.sources), ", "))

  iters, _, chains = size(c.value)
  value = Array(Float64, iters, length(nodenames), chains)
  for k in 1:chains
    for i in 1:iters
      relist!(m, vec(c.value[i,idx,k]), node.sources)
      v = isa(node.distr, Array) ? map(rand, node.distr) : rand(node.distr)
      value[i,:,k] = link(node, v, false)
    end
  end

  Chains(value, start=start(c.range), thin=step(c.range), names = nodenames)
end
