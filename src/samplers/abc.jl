#################### Approximate Bayesian Computing ####################

#################### Types and Constructors ####################

type ABCTune
  Told::Matrix{Float64}
  epsilon::Vector{Float64}
  terminals::Vector{Symbol}
end

ABCTune() = ABCTune(Matrix{Float64}(), Vector{Float64}(), Vector{Symbol}())

function ABC(params::Vector{Symbol}, sigma::Real, summarize::Vector{Function},
             target::Real; rho::Function = (x, y) -> sqrt(sumabs2(x - y)),
             s::Integer = 1, n::Integer = 50, kernel::Symbol = :uniform)

  kernel in [:uniform, :epanechnikov, :gaussian] ||
    throw(ArgumentError(
      "kernel must be one of :uniform, :epanechnikov, or :gaussian"
    ))

  samplerfx = function(model::Model, block::Integer)
    tunepar = tune(model, block)

    ## previous value
    x = unlist(model, block, true)
    logfx = mapreduce(p -> logpdf(model[p], true), +, params)
    

    T0 = NaN
    Told = NaN
    Tnew = NaN

    ## Initialize tune
    if model.iter == 1
      terminals = keys(model, :output)
      nt = length(terminals)

      ## data
      dvec = map(node -> unlist(model, [node]), terminals)

      ## summary stat of data
      T0vec = Array(Array{Float64, 1}, nt)
      for d in 1:nt
        T0vec[d] = vcat(map(f -> f(dvec[d]), summarize)...)
      end
      T0 = vcat(T0vec...)

      m = length(T0)

      Told = zeros(m, s)
      eps0 = zeros(s)
      for i in 1:s
        Toldvec = Array(Array{Float64, 1}, nt)
        for d in 1:nt
          x0 = unlist(model[terminals[d]], rand(model[terminals[d]]))
          Toldvec[d] = vcat(map(f -> f(x0), summarize)...)
        end
        Told[:, i] = vcat(Toldvec...)
        eps0[i] = rho(Told[:, i], T0)
      end
      tunepar.Told = Told
      tunepar.epsilon = eps0
      tunepar.terminals = terminals
    else
      nt = length(tunepar.terminals)
      ## data
      dvec = map(node -> unlist(model, [node]), tunepar.terminals)

      ## summary stat of data
      T0vec = Array(Array{Float64, 1}, nt)
      for d in 1:nt
        T0vec[d] = vcat(map(f -> f(dvec[d]), summarize)...)
      end
      T0 = vcat(T0vec...)
      Told = tunepar.Told
    end

    nt = length(tunepar.terminals)

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
      for i in 1:s
        ## new data
        Tnewvec = Array(Array{Float64, 1}, nt)
        for d in 1:nt
          x0 = unlist(model[tunepar.terminals[d]], rand(model[tunepar.terminals[d]]))
          Tnewvec[d] = vcat(map(f -> f(x0), summarize)...)
        end

        ## summary stat of new data
        Tnew[:, i] = vcat(Tnewvec...)

        ## Monotonically decrease epsilon to target
        epsilon[i] = max(target, min(rho(Tnew[:, i], T0), tunepar.epsilon[i]))

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
        tunepar.Told = Tnew
        tunepar.epsilon = epsilon
        break
      end
    end
    relist(model, x, block, true)
  end
  Sampler(params, samplerfx, ABCTune())
end
