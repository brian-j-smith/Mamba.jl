#################### MCMC Simulation Engine ####################

function mcmc{T<:String,U<:String}(model::MCMCModel, inputs::Dict{T},
           inits::Vector{Dict{U,Any}}, iter::Integer; burnin::Integer=0,
           thin::Integer=1, chains::Integer=1)

  iter > burnin || error("iter <= burnin")

  m = deepcopy(model)

  setinputs!(m, inputs)
  m.burnin = burnin
  tune0 = blocktune(m)

  monitorkeys = keys(m, true)

  sims = MCMCChain(labels(m, monitorkeys), div(iter - burnin - 1, thin) + 1,
                   start=burnin + thin, thin=thin, chains=chains, model=m)

  for k in 1:chains
    setinits!(m, inits[k])
    settune!(m, tune0)
    m.chain = k

    print("\nSAMPLING FROM CHAIN $(k)/$(chains)\n")
    pct = 0
    i = 1
    for t in 1:iter
      m.iter = t

      if floor(100 * t / iter) >= pct
        print(string("Iteration: ", lpad(t, length(string(iter)), ' '),
          "/$(iter) [", lpad(pct, 3, ' '), "%] @ $(strftime(time()))\n"))
        pct += 10
      end

      simulate!(m)

      if t > burnin && (t - burnin - 1) % thin == 0
        sims.data[i,:,k] = unlist(m, monitorkeys)
        i += 1
      end
    end
  end
  print("\n")

  sims
end
