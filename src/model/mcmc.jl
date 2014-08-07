#################### MCMC Simulation Engine ####################


function mcmc(c::Chains, iters::Integer)

  ismodelbased(c) || error("mcmc restart requires Chains from a Model fit")

  burnin = last(c.range) - c.model.iter
  thin = step(c.range)
  burnin + thin > 0 ||
    error("Chains have been subsetted to exclude the last iteration")

  mm = deepcopy(c.model)
  c2 = mcmc_master!(mm, mm.states[c.chains], iters, burnin, thin, size(c, 3))
  c2.model.iter += c.model.iter
  if c2.names != c.names
    c2 = c2[:,c.names,:]
  end

  Chains(cat(1, c.value, c2.value), start=first(c.range), thin=thin,
         names=c.names, model=c2.model)
end


function mcmc(m::Model, inputs::Dict{Symbol}, inits::Vector{Dict{Symbol,Any}},
              iters::Integer; burnin::Integer=0, thin::Integer=1,
              chains::Integer=1)

  iters > burnin || error("iters <= burnin")
  length(inits) >= chains || error("fewer initial values than chains")

  mm = deepcopy(m)
  setinputs!(mm, inputs)
  mm.states = Array(Vector{VariateType}, chains)
  mm.burnin = burnin

  mcmc_master!(mm, inits, iters, burnin, thin, chains)
end


function mcmc_master!(m::Model, inits, iters::Integer, burnin::Integer,
                      thin::Integer, chains::Integer)

  print(string("MCMC Simulation of $iters Iterations x $chains Chain",
               ifelse(chains > 1, "s", ""), "...\n\n"))

  lsts = [Any[m, inits[k], iters, burnin, thin, k] for k in 1:chains]
  sims = pmap(mcmc_worker!, lsts)

  m = sims[1].model
  m.states = map(k -> sims[k].model.states[k], 1:chains)
  Chains(cat(3, map(k -> sims[k].value, 1:chains)...),
         start=burnin+thin, thin=thin, names=sims[1].names, model=m)
end


function mcmc_worker!(args::Vector)

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
    if t > burnin && (t - burnin) % thin == 0
      sim.value[i,:,1] = unlist(model, true)
      i += 1
    end
    next!(meter)
  end
  model.states[chain] = unlist(model)
  println()

  sim
end
