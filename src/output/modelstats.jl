#################### Model-Based Posterior Statistics ####################

function dic(mc::ModelChains)
  nodekeys = keys(mc.model, :output)

  Dhat = -2.0 * logpdf(mc, mean, nodekeys)
  D = -2.0 * logpdf(mc, nodekeys).value
  p = [mean(D) - Dhat, 0.5 * var(D)]

  ChainSummary([Dhat + 2.0 * p  p], ["pD", "pV"],
               ["DIC", "Effective Parameters"], header(mc))
end


function logpdf(mc::ModelChains, f::Function, nodekeys::Vector{Symbol})
  m = mc.model

  relistkeys, updatekeys = getsimkeys(mc, nodekeys)
  relistkeys = union(relistkeys, intersect(nodekeys, keys(m, :block)))
  inds = names2inds(mc, relistkeys)

  m[relistkeys] = relist(m, map(i -> f(mc.value[:,i,:]), inds), relistkeys)
  update!(m, updatekeys)
  mapreduce(key -> logpdf(m[key]), +, nodekeys)
end


logpdf(mc::ModelChains, nodekey::Symbol) = logpdf(mc, [nodekey])

function logpdf(mc::ModelChains, nodekeys::Vector{Symbol})
  m = mc.model

  relistkeys, updatekeys = getsimkeys(mc, nodekeys)
  relistkeys = union(relistkeys, intersect(nodekeys, keys(m, :block)))
  inds = names2inds(mc, relistkeys)

  c = Chains(size(mc, 1), 1, chains=size(mc, 3), start=first(mc.range),
             thin=step(mc.range), names=["logpdf"])

  iters, _, chains = size(c.value)
  frame = ChainProgressFrame(
    "MCMC Processing of $iters Iterations x $chains Chain" * "s"^(chains > 1),
    true
  )
  for k in 1:chains
    meter = ChainProgress(frame, k, iters)
    for i in 1:iters
      m[relistkeys] = relist(m, mc.value[i,inds,k], relistkeys)
      update!(m, updatekeys)
      c.value[i,1,k] = mapreduce(key -> logpdf(m[key]), +, nodekeys)
      next!(meter)
    end
    println()
  end

  ModelChains(c, m)
end


predict(mc::ModelChains, nodekey::Symbol) = predict(mc, [nodekey])

function predict(mc::ModelChains, nodekeys::Vector{Symbol})
  m = mc.model

  outputs = keys(m, :output)
  all(key -> key in outputs, nodekeys) ||
    throw(ArgumentError(string(
      "nodekeys are not all observed Stochastic nodess : ",
      join(map(string, outputs), ", ")
    )))

  nodenames = names(m, nodekeys)
  relistkeys, updatekeys = getsimkeys(mc, nodekeys)
  inds = names2inds(mc, relistkeys)

  c = Chains(size(mc, 1), length(nodenames), chains=size(mc, 3),
             start=first(mc.range), thin=step(mc.range), names=nodenames)

  iters, _, chains = size(c.value)
  for k in 1:chains
    for i in 1:iters
      m[relistkeys] = relist(m, mc.value[i,inds,k], relistkeys)
      update!(m, updatekeys)
      f = key -> unlist(m[key], rand(m[key]))
      c.value[i,:,k] = vcat(map(f, nodekeys)...)
    end
  end

  ModelChains(c, m)
end


#################### Auxiliary Functions ####################

function getsimkeys(mc::ModelChains, nodekeys::Vector{Symbol})
  relistkeys = Symbol[]
  updatekeys = Symbol[]
  m = mc.model
  terminal = union(keys(m, :stochastic), keys(mc, :dependent))
  g = graph(m)
  for v in vertices(g)
    if isa(m[v.key], AbstractDependent)
      if any(key -> key in nodekeys, gettargets(v, g, m, terminal))
        v.key in terminal ?
          push!(relistkeys, v.key) :
          push!(updatekeys, v.key)
      end
    end
  end
  if !isempty(relistkeys) append!(updatekeys, nodekeys) end
  relistkeys, intersect(keys(m, :dependent), updatekeys)
end
