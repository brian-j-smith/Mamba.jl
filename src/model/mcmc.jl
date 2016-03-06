#################### MCMC Simulation Engine ####################

function mcmc(mc::ModelChains, iters::Integer; verbose::Bool=true)
  thin = step(mc)
  last(mc) == div(mc.model.iter, thin) * thin ||
    throw(ArgumentError("chain is missing its last iteration"))

  mm = deepcopy(mc.model)
  mc2 = mcmc_master!(mm, mm.iter + (1:iters), last(mc), thin, mc.chains,
                     verbose)
  if mc2.names != mc.names
    mc2 = mc2[:, mc.names, :]
  end

  ModelChains(vcat(mc, mc2), mc2.model)
end


function mcmc(m::Model, inputs::Dict{Symbol}, inits::Vector{Dict{Symbol, Any}},
              iters::Integer; burnin::Integer=0, thin::Integer=1,
              chains::Integer=1, verbose::Bool=true)
  iters > burnin ||
    throw(ArgumentError("burnin is greater than or equal to iters"))
  length(inits) >= chains ||
    throw(ArgumentError("fewer initial values than chains"))

  mm = deepcopy(m)
  setinputs!(mm, inputs)
  setinits!(mm, inits[1:chains])
  mm.burnin = burnin

  mcmc_master!(mm, 1:iters, burnin, thin, 1:chains, verbose)
end


function mcmc_master!(m::Model, window::UnitRange{Int}, burnin::Integer,
                      thin::Integer, chains::AbstractArray{Int}, verbose::Bool)
  states = m.states
  m.states = ModelState[]

  N = length(window)
  K = length(chains)

  frame = ChainProgressFrame(
    "MCMC Simulation of $N Iterations x $K Chain" * "s"^(K > 1), verbose
  )

  lsts = [
    Any[m, states[k], window, burnin, thin, ChainProgress(frame, k, N)]
    for k in chains
  ]
  results = pmap2(mcmc_worker!, lsts)

  sims  = Chains[results[k][1] for k in 1:K]
  model = results[1][2]
  model.states = ModelState[results[k][3] for k in sortperm(chains)]

  ModelChains(cat(3, sims...), model)
end


function mcmc_worker!(args::Vector)
  m, state, window, burnin, thin, meter = args

  m.iter = first(window) - 1
  relist!(m, state.value)
  settune!(m, state.tune)

  pnames = names(m, true)
  sim = Chains(last(window), length(pnames), start=burnin + thin, thin=thin,
               names=pnames)

  reset!(meter)
  for i in window
    sample!(m)
    if i > burnin && (i - burnin) % thin == 0
      sim[i, :, 1] = unlist(m, true)
    end
    next!(meter)
  end

  (sim, m, ModelState(unlist(m), gettune(m)))
end
