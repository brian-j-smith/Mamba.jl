#################### Approximate Bayesian Computation ####################

#################### Types and Constructors ####################

type ABCTune
  datakeys::Vector{Symbol}
  Tsim::Vector{Vector{Float64}}
  epsilon::Vector{Float64}

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
                      nsim::Integer=1, args...)

  params = asvec(params)
  kernelpdf = (epsilon, d) -> pdf(kernel(0.0, epsilon), d)
  local obsdata, simdata, summarizenodes

  samplerfx = function(model::Model, block::Integer)
    tune = Mamba.tune(model, block)

    ## current parameter and density values
    theta0 = unlist(model, block, true)
    logprior0 = mapreduce(key -> logpdf(model[key], true), +, params)
    pi_epsilon0 = 0.0

    ## initialize tuning parameters
    local Tobs
    if model.iter == 1
      ## data node symbols
      targets = keys(model, :target, params)
      stochastics = keys(model, :stochastic)
      tune.datakeys = intersect(setdiff(targets, params), stochastics)

      ## local functions to get and summarize data nodes
      obsdata = key -> unlist(model, [key])
      simdata = key -> unlist(model[key], rand(model[key]))
      summarizenodes = length(tune.datakeys) > 1 ?
        data -> vcat(map(key -> summary(data(key)), tune.datakeys)...) :
        data -> asvec(summary(data(tune.datakeys[1])))

      ## observed data summary statistics
      Tobs = summarizenodes(obsdata)

      tune.Tsim = Array(Vector{Float64}, nsim)
      tune.epsilon = Array(Float64, nsim)
      for i in 1:nsim
        ## simulated data summary statistics for current parameter values
        tune.Tsim[i] = summarizenodes(simdata)
        d = dist(tune.Tsim[i], Tobs; args...)

        ## starting tolerance
        tune.epsilon[i] = max(epsilon, d)

        ## kernel density evaluation
        pi_epsilon0 += kernelpdf(tune.epsilon[i], d)
      end
    else
      Tobs = summarizenodes(obsdata)
      for i in 1:nsim
        d = dist(tune.Tsim[i], Tobs; args...)
        pi_epsilon0 += kernelpdf(tune.epsilon[i], d)
      end
    end

    Tsim1 = similar(tune.Tsim)
    epsilon1 = similar(tune.epsilon)

    for k in 1:maxdraw
      ## candidate draw and prior density value
      theta1 = theta0 + scale .* rand(proposal(0.0, 1.0), length(theta0))
      relist!(model, theta1, block, true)
      logprior1 = mapreduce(key -> logpdf(model[key], true), +, params)

      ## tolerances and kernel density
      pi_epsilon1 = 0.0
      for i in 1:nsim
        ## simulated data summary statistics for candidate draw
        Tsim1[i] = summarizenodes(simdata)
        d = dist(Tsim1[i], Tobs; args...)

        ## monotonically decrease tolerance to target
        epsilon1[i] = max(epsilon, min(d, tune.epsilon[i]))

        ## kernel density evaluation
        pi_epsilon1 += kernelpdf(epsilon1[i], d)
      end

      ## accept/reject the candidate draw
      if rand() < pi_epsilon1 / pi_epsilon0 * exp(logprior1 - logprior0)
        theta0 = theta1
        tune.Tsim = Tsim1
        tune.epsilon = epsilon1
        break
      end
    end

    relist(model, theta0, block, true)
  end

  Sampler(params, samplerfx, ABCTune())
end
