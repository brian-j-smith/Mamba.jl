################## Binary Hamiltonian Monte Carlo  ##################

#################### Types ####################

type BHMCTune
  traveltime::Float64
  position::Vector{Float64}
  velocity::Vector{Float64}
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
    NaN,
    rand(Normal(0,1),length(x)),
    rand(Normal(0,1),length(x)),
    0,
    0
  )
  BHMCVariate(x, tune)
end

#################### Sampler Constructor ####################

function BHMC(params::Vector{Symbol}, traveltime::Real)
  Sampler(params,
  quote
    x = unlist(model,block)
    tunepar = tune(model, block)
    v = BHMCVariate(x,tunepar["sampler"])
    f = x -> logpdf!(model,x,block)
    bhmc!(v, tunepar["traveltime"], f)
    tunepar["sampler"] = v.tune
    relist(model,x,block)
  end,
  Dict("traveltime" => convert(Float64,traveltime), "sampler" => nothing)
  )
end

#################### Sampling Functions ####################

function bhmc!(v::BHMCVariate, traveltime::Float64, logf::Function)
  tune = v.tune
  flag = false
  nearzero = 10000*eps()
  j = 0
  totaltime = 0                       # records how much time the particle already moved

  d = size(v.value,1)                 # length of binary vector
  S = sign(tune.position)
  
  while true
    a = tune.velocity[:]
    b = tune.position[:]
    phi = atan2(b,a)

    walltime = -phi
    idx = find(x-> x>0,phi)                   
    walltime[idx] = pi - phi[idx]     # time to hit or cross wall

    # if there was a previous reflection (j>0)
    # and there is a potential reflection at the sample plane
    # make sure that a new reflection at j is not found because of numerical error
    if j > 0
      if abs(walltime[j]) < nearzero || abs(walltime[j] - 2*pi) < nearzero
        walltime[j] = Inf
      end
    end

    movetime, j = findmin(walltime)        # time till particle j hits or crosses wall
    if movetime == 0
      error("walking length zero!")
    elseif movetime == Inf
      movetime = pi
    end

    totaltime = totaltime + movetime
    if totaltime >= traveltime
      movetime = movetime - (totaltime-traveltime)
      flag = true
    else
      tune.wallhits += 1
    end

    # move the particle a time mt
    tune.position[:] = a*sin(movetime) + b*cos(movetime)
    tune.velocity[:] = a*cos(movetime) - b*sin(movetime)

    if flag
      break
    end

    tune.position[j] = 0

    S1 = (S + ones(d))/2; S1[j] = 0
    S2 = (S + ones(d))/2; S2[j] = 1

    v2_new = tune.velocity[j]^2 + sign(tune.velocity[j]) * 2 * (logf(S2)-logf(S1))
    if v2_new > 0
      tune.velocity[j] = sqrt(v2_new) * sign(tune.velocity[j])
      S[j] = -S[j]
      tune.wallcrosses += 1
    else
      tune.velocity[j] = -tune.velocity[j]
    end
  end
  v[:] = (sign(tune.position) + ones(d))/2 # convert from (-/+1) to (0/1)
end

