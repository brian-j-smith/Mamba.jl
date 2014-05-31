#################### MCMC Simulation Engine ####################

function mcmc(model::MCMCModel, inputs::Dict{Symbol},
           inits::Vector{Dict{Symbol,Any}}, iters::Integer; burnin::Integer=0,
           thin::Integer=1, chains::Integer=1)

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

    sims[k] = Array(VariateType, div(iters - burnin - 1, thin) + 1,
                    length(unlist(m, true)))

    print("\nSAMPLING FROM CHAIN $(k)/$(chains)\n")
    pct = 0
    i = 1
    for t in 1:iters
      m.iter = t

      if floor(100 * t / iters) >= pct
        print(string("Iteration: ", lpad(t, length(string(iters)), ' '),
          "/$(iters) [", lpad(pct, 3, ' '), "%] @ $(strftime(time()))\n"))
        pct += 10
      end

      simulate!(m)

      if t > burnin && (t - burnin - 1) % thin == 0
        sims[k][i,:] = unlist(m, true)
        i += 1
      end
    end
  end
  print("\n")

  MCMCChains(cat(3, sims...), names(m, true), start=burnin+thin, thin=thin,
             model=m)
end
