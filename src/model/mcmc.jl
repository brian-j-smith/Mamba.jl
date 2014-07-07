#################### MCMC Simulation Engine ####################

function mcmc(model::Model, inputs::Dict{Symbol},
              inits::Vector{Dict{Symbol,Any}}, iters::Integer;
              burnin::Integer=0, thin::Integer=1, chains::Integer=1)

  iters > burnin || error("iters <= burnin")

  m = deepcopy(model)
  setinputs!(m, inputs)

  sims = Array(Chains, chains)

  print(string("MCMC Simulation of $iters Iterations x $chains Chain",
               ifelse(chains > 1, "s", ""), "...\n\n"))

  np = nprocs()
  i = 1
  nextidx() = (idx = i; i += 1; idx)
  @sync begin
    for j in 1:np
      if j != myid() || np == 1
        @async begin
          while true
            k = nextidx()
            if k > chains
              break
            end
            sims[k] = remotecall_fetch(j, mcmc_sub!, m, inits[k], iters, burnin,
                                       thin, k)
          end
        end
      end
    end
  end

  cat(sims)
end

function mcmc_sub!(model::Model, inits::Dict{Symbol,Any}, iters::Integer,
                   burnin::Integer, thin::Integer, chain::Integer)

  setinits!(model, inits)
  model.burnin = burnin
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
  println()

  sim
end
