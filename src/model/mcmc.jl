#################### MCMC Simulation Engine ####################


function mcmc(mc::ModelChains, iters::Integer; verbose::Bool=true)
  burnin = last(mc.range) - mc.model.iter
  thin = step(mc.range)
  burnin + thin > 0 ||
    error("Chains have been subsetted to exclude the last iteration")

  mm = deepcopy(mc.model)
  mc2 = mcmc_master!(mm, mm.states[mc.chains], iters, burnin, thin, size(mc,3),
                    verbose)
  mc2.model.iter += mc.model.iter
  if mc2.names != mc.names
    mc2 = mc2[:,mc.names,:]
  end

  ModelChains(Chains(cat(1, mc.value, mc2.value), start=first(mc.range),
                     thin=thin, names=mc.names), mc2.model)
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
  ModelChains(Chains(cat(3, map(k -> sims[k].value, 1:chains)...),
                     start=burnin+thin, thin=thin, names=sims[1].names), m)
end


function mcmc_worker!(args::Vector)

  model, inits, iters, burnin, thin, chain, meter = args

  setinits!(model, inits)
  model.chain = chain

  pnames = names(model, true)
  sim = ModelChains(Chains(iters, length(pnames), start=burnin+thin, thin=thin,
                           names=pnames), model)

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
