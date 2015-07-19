#################### MCMC Simulation Engine ####################


function mcmc(c::Chains, iters::Integer; verbose::Bool=true)

  ismodelbased(c) || error("mcmc restart requires Chains from a Model fit")

  burnin = last(c.range) - c.model.iter
  thin = step(c.range)
  burnin + thin > 0 ||
    error("Chains have been subsetted to exclude the last iteration")

  mm = deepcopy(c.model)
  c2 = mcmc_master!(mm, mm.states[c.chains], iters, burnin, thin, size(c, 3),
                    verbose)
  c2.model.iter += c.model.iter
  if c2.names != c.names
    c2 = c2[:,c.names,:]
  end

  Chains(cat(1, c.value, c2.value), start=first(c.range), thin=thin,
         names=c.names, model=c2.model)
end


function mcmc(m::Model, inputs::Dict{Symbol}, inits::Vector{Dict{Symbol,Any}},
              iters::Integer; burnin::Integer=0, thin::Integer=1,
              chains::Integer=1, verbose::Bool=true)

  iters > burnin || error("iters <= burnin")
  length(inits) >= chains || error("fewer initial values than chains")

  mm = deepcopy(m)
  setinputs!(mm, inputs)
  mm.states = Array(Vector{Float64}, chains)
  mm.burnin = burnin

  mcmc_master!(mm, inits, iters, burnin, thin, chains, verbose)
end


function mcmc_master!(m::Model, inits, iters::Integer, burnin::Integer,
                      thin::Integer, chains::Integer, verbose::Bool)

  frame = ChainProgressFrame(
    "MCMC Simulation of $iters Iterations x $chains Chain" * "s"^(chains > 1),
    verbose
  )

  lsts = [
    Any[m, inits[k], iters, burnin, thin, k, ChainProgress(frame, k, iters)]
    for k in 1:chains
  ]
  sims = pmap(mcmc_worker!, lsts)

  m = sims[1].model
  m.states = map(k -> sims[k].model.states[k], 1:chains)
  Chains(cat(3, map(k -> sims[k].value, 1:chains)...),
         start=burnin+thin, thin=thin, names=sims[1].names, model=m)
end


function mcmc_worker!(args::Vector)

  model, inits, iters, burnin, thin, chain, meter = args

  setinits!(model, inits)
  model.chain = chain

  pnames = names(model, true)
  sim = Chains(iters, length(pnames), start=burnin+thin, thin=thin,
               names=pnames, model=model)

  reset!(meter)
  for i in 1:iters
    simulate!(model)
    if i > burnin && (i - burnin) % thin == 0
      sim[i,:,1] = unlist(model, true)
    end
    next!(meter)
  end
  model.states[chain] = unlist(model)
  println()

  sim
end
