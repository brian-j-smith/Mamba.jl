#################### Approximate Bayesian Computation ####################

#################### Types and Constructors ####################

type ABCTune
  datakeys::Vector{Symbol}
  Tsim::Vector{Vector{Float64}}
  epsilon::Vector{Float64}
  epsilonprime::Vector{Float64}

  function ABCTune()
    new(
      Vector{Symbol}(),
      Vector{Vector{Float64}}(),
      Vector{Float64}()
    )
  end
end


#################### Sampler Constructor ####################

function ABC{T<:Real}(params::ElementOrVector{Symbol},
                      scale::ElementOrVector{T}, summary::Function,
                      epsilon::Real; kernel::KernelDensityType=SymUniform,
                      dist::Function=(Tsim, Tobs) -> sqrt(sumabs2(Tsim - Tobs)),
                      proposal::SymDistributionType=Normal, maxdraw::Integer=1,
                      nsim::Integer=1, decay::Real=1.0, randeps::Bool=false,
                      args...)
  0 <= decay <= 1 || throw(ArgumentError("decay is not in [0, 1]"))

  params = asvec(params)
  kernelpdf = (epsilon, d) -> pdf(kernel(0.0, epsilon), d)
  local obsdata, simdata, summarizenodes

  samplerfx = function(model::Model, block::Integer)
    tune = gettune(model, block)

    ## current parameter and density values
    theta0 = unlist(model, block, true)
    logprior0 = logpdf(model, params, true)
    pi_epsilon0 = 0.0
    pi_error0 = 1.0

    ## initialize tuning parameters
    local Tobs
    if model.iter == 1
      ## data node symbols
      targets = keys(model, :target, params)
      stochastics = keys(model, :stochastic)
      tune.datakeys = intersect(setdiff(targets, params), stochastics)

      ## local functions to get and summarize data nodes
      obsdata = key -> unlist(model[key])
      simdata = key -> unlist(model[key], rand(model[key]))
      summarizenodes = length(tune.datakeys) > 1 ?
        data -> vcat(map(key -> summary(data(key)), tune.datakeys)...) :
        data -> asvec(summary(data(tune.datakeys[1])))

      ## observed data summary statistics
      Tobs = summarizenodes(obsdata)

      tune.Tsim = Array{Vector{Float64}}(nsim)
      tune.epsilon = Array{Float64}(nsim)
      tune.epsilonprime = Array{Float64}(nsim)
      for i in 1:nsim
        ## simulated data summary statistics for current parameter values
        tune.Tsim[i] = summarizenodes(simdata)
        d = dist(tune.Tsim[i], Tobs; args...)

        ## starting tolerance
        tune.epsilon[i] = decay > 0 ? max(epsilon, d) : epsilon
        if randeps
          dexp = Exponential(tune.epsilon[i])
          tune.epsilonprime[i] = rand(dexp)
          pi_error0 = pdf(dexp, tune.epsilonprime[i])
        else
          tune.epsilonprime[i] = tune.epsilon[i]
        end

        ## kernel density evaluation
        pi_epsilon0 += kernelpdf(tune.epsilonprime[i], d) * pi_error0
      end
    else
      Tobs = summarizenodes(obsdata)
      for i in 1:nsim
        d = dist(tune.Tsim[i], Tobs; args...)
        if randeps
          dexp = Exponential(tune.epsilon[i])
          pi_error0 = pdf(dexp, tune.epsilonprime[i])
        end
        pi_epsilon0 += kernelpdf(tune.epsilonprime[i], d) * pi_error0
      end
    end

    Tsim1 = similar(tune.Tsim)
    epsilon1 = similar(tune.epsilon)
    epsilonprime1 = similar(tune.epsilonprime)

    for k in 1:maxdraw
      ## candidate draw and prior density value
      theta1 = theta0 + scale .* rand(proposal(0.0, 1.0), length(theta0))
      relist!(model, theta1, block, true)
      logprior1 = logpdf(model, params, true)

      ## tolerances and kernel density
      pi_epsilon1 = 0.0
      pi_error1 = 1.0
      for i in 1:nsim
        ## simulated data summary statistics for candidate draw
        Tsim1[i] = summarizenodes(simdata)
        d = dist(Tsim1[i], Tobs; args...)

        ## monotonically decrease tolerance to target
        epsilon1[i] = (1 - decay) * tune.epsilon[i] +
                      decay * max(epsilon, min(d, tune.epsilon[i]))
        if randeps
          dexp = Exponential(epsilon1[i])
          epsilonprime1[i] = rand(dexp)
          pi_error1 = pdf(dexp, epsilonprime1[i])
        else
          epsilonprime1[i] = epsilon1[i]
        end

        ## kernel density evaluation
        pi_epsilon1 += kernelpdf(epsilonprime1[i], d) * pi_error1
      end

      ## accept/reject the candidate draw
      if rand() < pi_epsilon1 / pi_epsilon0 * exp(logprior1 - logprior0)
        theta0 = theta1
        tune.Tsim = Tsim1
        tune.epsilon = epsilon1
        tune.epsilonprime = epsilonprime1
        break
      end
    end

    relist(model, theta0, block, true)
  end

  Sampler(params, samplerfx, ABCTune())
end
