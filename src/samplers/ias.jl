##################### Individual Adaptation ####################

#################### Types and Constructors ####################

type IASTune <: SamplerTune
  logf::Nullable{Function}
  A::Vector{Float64}
  D::Vector{Float64}

  epsilon::Float64
  lambda::Float64
  tau::Float64

  iter::Integer

  IASTune() = new()

  function IASTune(x::Vector, logf::Nullable{Function}; 
          A::Vector{Float64} = [1/length(x) for j in 1:length(x)],
          D::Vector{Float64} = [1/length(x) for j in 1:length(x)],
          epsilon::Float64 = 0.01 / length(x), lambda::Float64 = 0.55,
          tau::Float64 = 0.45)

    new(logf, A, D, epsilon, lambda, tau, 0)
  end
end

IASTune(x::Vector; args...) =
   IASTune(x, Nullable{Function}(); args...)

IASTune(x::Vector, logf::Function; args...) =
   IASTune(x, Nullable{Function}(logf); args...)


typealias IASVariate SamplerVariate{IASTune}

function validate(v::IASVariate)
  n = length(v)
  0.0 < v.tune.epsilon < 0.5 || throw(ArgumentError("epsilon must be 0.0 < epsilon < 0.5"))
  0.5 < v.tune.lambda <= 1.0 || throw(ArgumentError("lambda must be 0.5 < lambda <= 1.0"))
  0.0 < v.tune.tau < 1.0 || throw(ArgumentError("tau must be 0.0 < tau < 1.0"))
  validatebinary(v)
end


#################### Sampler Constructor ####################

function IAS(params::ElementOrVector{Symbol}; args...)
  samplerfx = function(model::Model, block::Integer)
    block = SamplingBlock(model, block)
    v = SamplerVariate(block; args...)
    sample!(v, x -> logpdf!(block, x))
    relist(block, v)
  end
  Sampler(params, samplerfx, IASTune())
end


#################### Sampling Functions ####################

sample!(v::IASVariate) = sample!(v, v.tune.logf)

function sample!(v::IASVariate, logf::Function)
  v.tune.iter += 1
  A_new = similar(v.tune.A)
  D_new = similar(v.tune.D)

  ## proposal
  x = v[:]
  p = length(x)
  log_q_num = 0.0
  log_q_den = 0.0 
  for j in 1:p
    if v.value[j] == 0
      if rand() < v.tune.A[j]
        x[j] = 1
        log_q_den += log(v.tune.A[j])
        log_q_num += log(v.tune.D[j])
      else
        log_q_den += log(1 - v.tune.A[j])
        log_q_num += log(1 - v.tune.A[j])
      end
    else
      if rand() < v.tune.D[j]
        x[j] = 0
        log_q_den += log(v.tune.D[j])
        log_q_num += log(v.tune.A[j])
      else
        log_q_den += log(1 - v.tune.D[j])
        log_q_num += log(1 - v.tune.D[j])
      end
    end
  end

  ## M-H acceptance probability
  alpha_p = exp(logf(x) - logf(v.value))
  alpha = min(1, exp( alpha_p + (log_q_num - log_q_den)))

  ## Adapt tuning parameters
  for j in 1:p
    added = 0.0
    deleted = 0.0

    if x[j] != v.value[j]
      if v.value[j] == 0
        added = 1.0
      else
        deleted=1.0
      end
    end

    ## new A[j]
    C = log((v.tune.A[j] - v.tune.epsilon) / (1 - v.tune.A[j] - v.tune.epsilon)) + 
        v.tune.iter ^ (-v.tune.lambda) * added * (alpha - v.tune.tau)
    A_new[j] = (exp(C) - v.tune.epsilon * exp(C) + v.tune.epsilon) / (1 + exp(C))
    
    ## new D[j]
    C = log((v.tune.D[j] - v.tune.epsilon) / (1 - v.tune.D[j] - v.tune.epsilon)) + 
        v.tune.iter ^ (-v.tune.lambda) * deleted * (alpha - v.tune.tau)
    D_new[j] = (exp(C) - v.tune.epsilon * exp(C) + v.tune.epsilon) / (1 + exp(C))
  end

  ## Accept or reject proposal
  if rand() < alpha
    v[:] = x 
  end

  v.tune.A[:] = A_new
  v.tune.D[:] = D_new
  v
end
