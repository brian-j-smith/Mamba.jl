#################### MCMC Simulation Engine ####################

function mcmc(model::Model, inputs::Dict{Symbol},
              inits::Vector{Dict{Symbol,Any}}, iters::Integer;
              burnin::Integer=0, thin::Integer=1, chains::Integer=1)

  iters > burnin || error("iters <= burnin")

  m = deepcopy(model)

  setinputs!(m, inputs)
  m.burnin = burnin
  tune0 = tune(m)

  sims = Array(Matrix{VariateType}, chains)

  for k in 1:chains
    setinits!(m, inits[k])
    settune!(m, tune0)
    m.chain = k

    sims[k] = Array(VariateType, length(burnin+1:thin:iters),
                    length(unlist(m, true)))

    i = 1
    meter = ChainProgress(k, iters)
    for t in 1:iters
      simulate!(m)
      if t > burnin && (t - burnin - 1) % thin == 0
        sims[k][i,:] = unlist(m, true)
        i += 1
      end
      next!(meter)
    end
    println()
  end

  Chains(cat(3, sims...), start=burnin+thin, thin=thin, names=names(m, true),
         model=m)
end
