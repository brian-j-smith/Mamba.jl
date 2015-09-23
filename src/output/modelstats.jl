#################### Model-Based Posterior Statistics ####################

function dic(mc::ModelChains)
  m = mc.model
  nkeys = keys(m, :output)
  idx = indexin(names(m, keys(m, :block)), mc.names)
  in(0, idx) && error("dic requires all sampled nodes to be monitored")

  xbar = map(i -> mean(mc.value[:,i,:]), idx)
  relist!(m, xbar)
  Dhat = -2.0 * mapreduce(key -> logpdf(m[key]), +, nkeys)
  D = -2.0 * logpdf(mc, nkeys)
  p = [mean(D) - Dhat, 0.5 * var(D)]

  ChainSummary([Dhat + 2.0 * p  p], ["pD", "pV"],
               ["DIC", "Effective Parameters"], header(mc))
end

function logpdf(mc::ModelChains, nkeys::Vector{Symbol})
  m = mc.model
  idx = indexin(names(m, keys(m, :block)), mc.names)
  in(0, idx) && error("logpdf requires all sampled nodes to be monitored")

  iters, p, chains = size(mc.value)
  values = Array(Float64, iters, 1, chains)
  frame = ChainProgressFrame(
    "MCMC Processing of $iters Iterations x $chains Chain" * "s"^(chains > 1),
    true
  )
  for k in 1:chains
    meter = ChainProgress(frame, k, iters)
    for i in 1:iters
      relist!(m, vec(mc.value[i,idx,k]))
      values[i,1,k] = mapreduce(key -> logpdf(m[key]), +, nkeys)
      next!(meter)
    end
    println()
  end

  values
end

function predict(mc::ModelChains, key::Symbol)
  m = mc.model
  node = m[key]

  outputs = keys(m, :output)
  in(key, outputs) ||
    error("predict is only defined for observed Stochastic nodes: ",
          join(map(string, outputs), ", "))

  nodenames = names(m, [key])

  sources = intersect(node.sources, keys(m, :stochastic))
  idx = indexin(names(m, sources), mc.names)
  in(0, idx) && error("predict requires monitoring of nodes: ",
                      join(map(string, sources), ", "))

  iters, _, chains = size(mc.value)
  value = Array(Float64, iters, length(nodenames), chains)
  for k in 1:chains
    for i in 1:iters
      relist!(m, vec(mc.value[i,idx,k]), sources)

      if isa(node.distr, Distribution)
        x = rand(node.distr)
      elseif isa(node.distr, Array{UnivariateDistribution})
        x = map(rand, node.distr)
      elseif isa(node.distr, Array{MultivariateDistribution})
        x = Array(Float64, dims(node.distr))
        for sub in CartesianRange(size(node.distr))
          d = node.distr[sub]
          x[sub, 1:length(d)] = rand(d)
        end
      end

      value[i,:,k] = unlist(node, x)
    end
  end

  Chains(value, start=start(mc.range), thin=step(mc.range), names = nodenames)
end
