type ABCTune
  Told::Matrix{Float64}
  epsilon::Vector{Float64}
  terminals::Vector{Symbol}
end

function ABC(params::Vector{Symbol}, sigma::Real, summarize::Vector{Function},
             epsilon::Real; rho::Function = (x, y) -> sqrt(sum((x - y) .^ 2)),
             s::Int = 1, n::Int = 50, kernel::Symbol = :uniform)
  kernel in [:uniform, :epanechnikov, :gaussian] ||
    throw(ArgumentError("adapt must be one of :uniform, 
    :epanechnikov, or :gaussian"))

  Sampler(params,
    quote
      tunepar = tune(model, block)
      params = keys(model, :block, block)

      #previous value
      x = unlist(model, block, true)
      
      #logprior
      logprior = x -> begin
        value = 0.0
        relist!(model, x, block, true)
        for key in params
          value += logpdf(model[key], true)
        end
        return value
      end

      T0 = NaN
      Told = NaN
      Tnew = NaN

      if tunepar["sampler"] == nothing 
        terminals = keys(model, :output)

        #data
        d = unlist(model, terminals)

        #summary stat of data
        T0 = vcat(map(f -> f(d), tunepar["summarize"])...)

        m = length(T0)

        Told = zeros(m, tunepar["s"])
        eps0 = zeros(tunepar["s"])
        for i in 1:tunepar["s"]
          x0 =  vcat(map(d -> unlist(model[d], 
                  rand(model[d])), terminals)...)
          Told[:, i] = vcat(map(f -> f(x0), tunepar["summarize"])...)
          eps0[i] = tunepar["rho"](Told[:, i], T0)
        end
        tunepar["sampler"] = ABCTune(Told, eps0, terminals)
      else
        #data
        d = unlist(model, tunepar["sampler"].terminals)

        #summary stat of data
        T0 = vcat(map(f -> f(d), tunepar["summarize"])...)
        Told = tunepar["sampler"].Told
      end

      for k in 1:tunepar["n"]
        #proposal
        y = x + tunepar["sigma"] * randn(length(x))
        relist!(model, y, block, true)

        ratio_num = 0.0
        ratio_den = 0.0

        Tnew = zeros(length(T0), tunepar["s"])
        epsilon = zeros(tunepar["s"])

        for i in 1:tunepar["s"]
          #new data
          dnew = vcat(map(d -> unlist(model[d], 
                   rand(model[d])), tunepar["sampler"].terminals)...)

          #summary stat of new data
          Tnew[:, i] = vcat(map(f -> f(dnew), tunepar["summarize"])...)

          #epsilon
          epsilon[i] = max(tunepar["epsilon"], min(tunepar["rho"](Tnew[:, i], 
                         T0), tunepar["sampler"].epsilon[i]))

          #weighting density
          rho0 = 0
          rho1 = 0
          if tunepar["kernel"] == :uniform
            rho0 = pdf(Uniform(0, epsilon[i]), tunepar["rho"](Told[:, i], T0))
            rho1 = pdf(Uniform(0, epsilon[i]), tunepar["rho"](Tnew[:, i], T0))
          elseif tunepar["kernel"] == :epanechnikov
            rho0 = prod([pdf(Epanechnikov(Told[j, i], epsilon[i]), T0[i]) 
              for j in 1:length(T0)])
            rho1 = prod([pdf(Epanechnikov(Tnew[j, i], epsilon[i]), T0[i]) 
              for j in 1:length(T0)])
          elseif tunepar["kernel"] == :gaussian
            rho0 = pdf(MvNormal(collect(Told[:, i]), epsilon[i]), collect(T0))
            rho1 = pdf(MvNormal(collect(Tnew[:, i]), epsilon[i]), collect(T0))
          end

          ratio_num += rho1
          ratio_den += rho0
        end


        #Reject/Accept
        if rand() < ratio_num / ratio_den * exp(logprior(y) - logprior(x))
          x[:] = y
          tunepar["sampler"].Told = Tnew
          tunepar["sampler"].epsilon = epsilon
          break
        end
      end
      relist(model, x, block, true)
    end,
    Dict("summarize" => summarize, "rho" => rho, "n" => convert(Int64, n), 
         "s" => convert(Int64, s), "epsilon" => convert(Float64, epsilon), 
         "sigma" => convert(Float64, sigma), "kernel" => kernel, 
         "sampler" => nothing)
    )
end
