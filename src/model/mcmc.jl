#################### MCMC Simulation Engine ####################

function mcmc(mc::ModelChains, iters::Integer; verbose::Bool=true)
  last(mc) > mc.model.iter - step(mc) ||
    throw(ArgumentError("chain is missing its last iteration"))

  mm = deepcopy(mc.model)
  mc2 = mcmc_master!(mm, mm.states[mc.chains], mm.iter + (1:iters), last(mc),
                     step(mc), size(mc, 3), verbose)
  mc2.model.iter += mc.model.iter
  if mc2.names != mc.names
    mc2 = mc2[:,mc.names,:]
  end

  ModelChains(vcat(mc, mc2), mc2.model)
end


function mcmc(m::Model, inputs::Dict{Symbol}, inits::Vector{Dict{Symbol,Any}},
              iters::Integer; burnin::Integer=0, thin::Integer=1,
              chains::Integer=1, verbose::Bool=true)
  @everywhere using Mamba

  iters > burnin ||
    throw(ArgumentError("burnin is greater than or equal to iters"))
  length(inits) >= chains ||
    throw(ArgumentError("fewer initial values than chains"))

  mm = deepcopy(m)
  setinputs!(mm, inputs)
  mm.states = Array(Vector{Float64}, chains)
  mm.burnin = burnin

  mcmc_master!(mm, inits, 1:iters, burnin, thin, chains, verbose)
end


function mcmc_master!(m::Model, inits, window::UnitRange{Int}, burnin::Integer,
                      thin::Integer, chains::Integer, verbose::Bool)
  n = length(window)
  frame = ChainProgressFrame(
    "MCMC Simulation of $n Iterations x $chains Chain" * "s"^(chains > 1),
    verbose
  )
  lsts = [
    Any[m, inits[k], window, burnin, thin, k, ChainProgress(frame, k, n)]
    for k in 1:chains
  ]
  sims = mcmcmap(mcmc_worker!, lsts)
  model = sims[1].model
  model.states = map(k -> sims[k].model.states[k], 1:chains)

  ModelChains(cat(3, sims...), model)
end


function mcmc_worker!(args::Vector)
  model, inits, window, burnin, thin, chain, meter = args

  setinits!(model, inits)
  model.iter = first(window) - 1
  model.chain = chain

  pnames = names(model, true)
  sim = ModelChains(Chains(last(window), length(pnames), start=burnin+thin,
                           thin=thin, names=pnames), model)

  reset!(meter)
  for i in window
    simulate!(model)
    if i > burnin && (i - burnin) % thin == 0
      sim[i,:,1] = unlist(model, true)
    end
    next!(meter)
  end
  model.states[chain] = unlist(model)

  sim
end


#################### Auxiliary Functions ####################

## mcmcmap is a partial work-around for the pmap issue in julia 0.4.0 of worker
## node errors being blocked.  In single-processor mode, mcmcmap calls map
## instead to avoid the error handling issue.  In multi-processor model, pmap is
## called and will apply its error processing.  If and when the pmap issue is
## resolved in a future version of julia, calls to mcmcmap should be reverted to
## a call to pmap.

function mcmcmap(f::Function, lsts::AbstractArray)
  nprocs() > 1 ? pmap(f, lsts) : map(f, lsts)
end
