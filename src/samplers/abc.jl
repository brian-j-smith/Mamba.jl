type ABCTune
  Told::Matrix{Float64}
  epsilon::Vector{Float64}
  terminals::Vector{Symbol}
end

function ABC(params::Vector{Symbol}, sigma::Real, summarize::Vector{Function},
             epsilon::Real; rho::Function = (x, y) -> sqrt(sumabs2(x - y)),
             s::Integer = 1, n::Integer = 50, kernel::Symbol = :uniform)

  kernel in [:uniform, :epanechnikov, :gaussian] ||
    throw(ArgumentError(
      "kernel must be one of :uniform, :epanechnikov, or :gaussian"
    ))

  Sampler(params, (model::Model, block::Integer) ->
    begin
      tunepar = tune(model, block)
      params = keys(model, :block, block)

      ## previous value
      x = unlist(model, block, true)
      logfx = mapreduce(p -> logpdf(model[p], true), +, params)
      
      T0 = NaN
      Told = NaN
      Tnew = NaN

      ## User inputs
      sigma = tunepar["sigma"]
      summarize = tunepar["summarize"]
      target = tunepar["epsilon"]
      rho = tunepar["rho"]
      s = tunepar["s"]
      n = tunepar["n"]
      kernel = tunepar["kernel"]

      if tunepar["sampler"] == nothing 
        terminals = keys(model, :output)

        ## data
        d = unlist(model, terminals)

        ## summary stat of data
        T0 = vcat(map(f -> f(d), summarize)...)

        m = length(T0)

        Told = zeros(m, s)
        eps0 = zeros(s)
        for i in 1:s
          x0 =  vcat(map(d -> unlist(model[d], rand(model[d])), 
                         terminals)...)
          Told[:, i] = vcat(map(f -> f(x0), summarize)...)
          eps0[i] = rho(Told[:, i], T0)
        end
        tunepar["sampler"] = ABCTune(Told, eps0, terminals)
      else
        ## data
        d = unlist(model, tunepar["sampler"].terminals)

        ## summary stat of data
        T0 = vcat(map(f -> f(d), summarize)...)
        Told = tunepar["sampler"].Told
      end

      ## Attempt n times to get new proposal
      for k in 1:n
        ## proposal
        y = x + sigma * randn(length(x))
        relist!(model, y, block, true)
        logfy = mapreduce(p -> logpdf(model[p], true), +, params)

        ratio_num = 0.0
        ratio_den = 0.0

        Tnew = zeros(length(T0), s)
        epsilon = zeros(s)

        ## Simulate s data sets
        for i in 1:tunepar["s"]
          ## new data
          dnew = vcat(map(d -> unlist(model[d], rand(model[d])), 
                          tunepar["sampler"].terminals)...)

          ## summary stat of new data
          Tnew[:, i] = vcat(map(f -> f(dnew), summarize)...)

          ## Monotonically decrease epsilon to target
          epsilon[i] = max(target, min(rho(Tnew[:, i], T0), 
                                       tunepar["sampler"].epsilon[i]))
          ## weighting density
          rho0 = 0
          rho1 = 0
          if kernel == :uniform
            rho0 = pdf(Uniform(0, epsilon[i]), rho(Told[:, i], T0))
            rho1 = pdf(Uniform(0, epsilon[i]), rho(Tnew[:, i], T0))
          elseif kernel == :epanechnikov
            rho0 = prod([pdf(Epanechnikov(Told[j, i], epsilon[i]), T0[j]) 
              for j in 1:length(T0)])
            rho1 = prod([pdf(Epanechnikov(Tnew[j, i], epsilon[i]), T0[j]) 
              for j in 1:length(T0)])
          elseif kernel == :gaussian
            rho0 = pdf(MvNormal(collect(Told[:, i]), epsilon[i]), collect(T0))
            rho1 = pdf(MvNormal(collect(Tnew[:, i]), epsilon[i]), collect(T0))
          end

          ratio_num += rho1
          ratio_den += rho0
        end

        ## Reject/Accept
        if rand() < ratio_num / ratio_den * exp(logfy - logfx)
          x[:] = y
          tunepar["sampler"].Told = Tnew
          tunepar["sampler"].epsilon = epsilon
          break
        end
      end
      relist(model, x, block, true)
    end,
    Dict("sigma" => convert(Float64, sigma), "summarize" => summarize, 
         "epsilon" => convert(Float64, epsilon), "rho" => rho, 
         "s" => convert(Int64, s), "n" => convert(Int64, n), 
         "kernel" => kernel, "sampler" => nothing)
    )
end
