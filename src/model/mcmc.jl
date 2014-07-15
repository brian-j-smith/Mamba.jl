#################### MCMC Simulation Engine ####################


function mcmc(model::Model, inputs::Dict{Symbol},
              inits::Vector{Dict{Symbol,Any}}, iters::Integer;
              burnin::Integer=0, thin::Integer=1, chains::Integer=1)

  iters > burnin || error("iters <= burnin")

  m = deepcopy(model)
  setinputs!(m, inputs)
  m.state = Array(Vector{VariateType}, chains)
  m.burnin = burnin

  print(string("MCMC Simulation of $iters Iterations x $chains Chain",
               ifelse(chains > 1, "s", ""), "...\n\n"))

  lsts = [[m, inits[k], iters, burnin, thin, k] for k in 1:chains]
  sims = pmap(mcmc_sub!, lsts)

  m = sims[1].model
  m.state = map(k -> sims[k].model.state[k], 1:chains)
  Chains(cat(3, map(k -> sims[k].value, 1:chains)...),
         start=burnin+thin, thin=thin, names=sims[1].names, model=m)
end


function mcmc_sub!(args::Vector)

  model, inits, iters, burnin, thin, chain = args

  setinits!(model, inits)
  model.chain = chain

  pnames = names(model, true)
  sim = Chains(iters, length(pnames), start=burnin+thin, thin=thin,
               names=pnames, model=model)

  i = 1
  meter = ChainProgress(chain, iters)
  for t in 1:iters
    simulate!(model)
    if t > burnin && (t - burnin - 1) % thin == 0
      sim.value[i,:,1] = unlist(model, true)
      i += 1
    end
    next!(meter)
  end
  model.state[chain] = unlist(model)
  println()

  sim
end
