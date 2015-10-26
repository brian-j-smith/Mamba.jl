function ABC(params::Vector{Symbol}, summarize::Function, 
             rho::Function, epsilon::Vector{Float64}, sigma::Float64; n::Int=50)
  Sampler(params,
    quote
      tunepar = tune(model, block)
      params = keys(model, :block, block)

      #logprior
      logprior = x -> begin
        value = 0.0
        values = relist(model, x, params, true)
        for key in params
          value += logpdf(model[key], values[key], true)
        end
        return value
      end

      #previous value
      x = unlist(model, block, true)
      
      #data
      d = unlist(model, keys(model, :output))

      #summary stat of data
      T0 = tunepar["summarize"](d)

      for i in 1:tunepar["n"]
        #proposal
        y = x + tunepar["sigma"] * randn(length(x))
        relist!(model, y, block, true)

        #new data
        dnew = vcat(map(d -> unlist(model[d], rand(model[d])), keys(model, :output))...)

        #summary stat of new data
        Tnew = tunepar["summarize"](dnew)

        #Reject/Accept
        if all(tunepar["rho"](Tnew, T0) .< tunepar["epsilon"])
          if rand() < exp(logprior(y)/logprior(x))
            x[:] = y
            break
          end
        end
      end

      relist(model, x, block, true)
    end,
    Dict("summarize" => summarize, "rho" => rho, "n" => n,
         "epsilon" => epsilon, "sigma" => sigma)
    )
end
