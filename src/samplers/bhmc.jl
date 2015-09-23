################## Binary Hamiltonian Monte Carlo  ##################

#################### Types ####################

type BHMCTune
  T::Float64
  X::Vector{Float64}
  V::Vector{Float64}
  wallhits::Int
  wallcrosses::Int
end

type BHMCVariate <: VectorVariate
  value::Vector{Float64}
  tune::BHMCTune

  BHMCVariate(x::Vector{Float64}, tune::BHMCTune) = new(x, tune)
end

function BHMCVariate(x::Vector{Float64}, tune=nothing)
  tune = BHMCTune(
    0,
    rand(Normal(0,1),length(x)),
    rand(Normal(0,1),length(x)),
    0,
    0
  )
  BHMCVariate(x, tune)
end

#################### Sampler Constructor ####################

function BHMC(params::Vector{Symbol}, T::Real)
  Sampler(params,
  quote
    x = unlist(model,block)
    tunepar = tune(model, block)
    v = BHMCVariate(x,tunepar["sampler"])
    f = x -> logpdf!(model,x,block)
    bhmc!(v, tunepar["T"], f)
    tunepar["sampler"] = v.tune
    relist(model,x,block)
  end,
  Dict("T" => convert(Float64,T), "sampler" => nothing)
  )
end

#################### Sampling Functions ####################

function bhmc!(v::BHMCVariate, T::Float64, logf::Function)
  tune = v.tune
  flag = false
  nearzero = 10000*eps()
  j = 0
  tt = 0                       # records how much time the particle already moved

  d = size(v.value,1)          # length of binary vector
  tune.V = rand(Normal(0,1),d) # initialize velocity/momentum

  S = sign(tune.X)

  while true
    a = tune.V[:]
    b = tune.X[:]
    phi = atan2(b,a)

    wt1 = -phi
    idx = find(phi)       # which phi are positive
    wt1[idx] = pi - phi[idx]

    # if there was a previous reflection (j>0)
    # and there is a potential reflection at the sample plane
    # make sure that a new reflection at j is not found because of numerical error
    if j > 0
      tt1 = wt1[j]
      if abs(tt1) < nearzero || abs(tt1 - 2*pi) < nearzero
        wt1[j] = Inf
      end
    end
    mt, j = findmin(wt1)
    if mt == 0
      error("walking length zero!")
    elseif mt == Inf
      mt = pi
    end

    tt = tt + mt
    if tt >= T
      mt = mt - (tt-T)
      flag = true
    else
      tune.wallhits += 1
    end

    # move the particle a time mt
    tune.X[:] = a*sin(mt) + b*cos(mt)
    tune.V[:] = a*cos(mt) - b*sin(mt)

    if flag
      break
    end

    tune.X[j] = 0

    S1 = (S + ones(d))/2; S1[j] = 0
    S2 = (S + ones(d))/2; S2[j] = 1

    v2_new = tune.V[j]^2 + sign(tune.V[j]) * 2 * (logf(S1)-logf(S2))
    if v2_new > 0
      tune.V[j] = sqrt(v2_new) * sign(tune.V[j])
      S[j] = -S[j]
      tune.wallcrosses += 1
    else
      tune.V[j] = -tune.V[j]
    end
  end
  v[:] = (sign(tune.X) + ones(d))/2 # convert from (-/+1) to (0/1)
end

