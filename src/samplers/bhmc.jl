################## Binary Hamiltonian Monte Carlo  ##################

#################### Types ####################

mutable struct BHMCTune <: SamplerTune
  logf::Union{Function, Missing}
  traveltime::Float64
  position::Vector{Float64}
  velocity::Vector{Float64}
  wallhits::Int
  wallcrosses::Int

  BHMCTune() = new()

  function BHMCTune(x::Vector, traveltime::Real, logf::Union{Function, Missing})
    n = length(x)
    new(logf, traveltime, randn(n), randn(n), 0, 0)
  end
end

BHMCTune(x::Vector, traveltime::Real) =
  BHMCTune(x, traveltime, missing)

const BHMCVariate = SamplerVariate{BHMCTune}

validate(v::BHMCVariate) = validatebinary(v)


#################### Sampler Constructor ####################

function BHMC(params::ElementOrVector{Symbol}, traveltime::Real)
  samplerfx = function(model::Model, block::Integer)
    block = SamplingBlock(model, block)
    v = SamplerVariate(block, traveltime)
    sample!(v, x -> logpdf!(block, x))
    relist(block, v)
  end
  Sampler(params, samplerfx, BHMCTune())
end


#################### Sampling Functions ####################

sample!(v::BHMCVariate) = sample!(v, v.tune.logf)

function sample!(v::BHMCVariate, logf::Function)
  tune = v.tune
  flag = false
  nearzero = 1e4 * eps()
  j = 0
  totaltime = 0.0                     ## time the particle already moved

  n = length(v)                       ## length of binary vector
  S = sign.(tune.position)

  while true
    a = tune.velocity[:]
    b = tune.position[:]
    phi = atan.(b, a)

    ## time to hit or cross wall
    walltime = -phi
    idx = findall(x-> x > 0.0, phi)
    walltime[idx] .= pi .- phi[idx]

    ## if there was a previous reflection (j > 0) and there is a potential
    ## reflection at the sample plane make sure that a new reflection at j
    ## is not found because of numerical error
    if j > 0
      if abs(walltime[j]) < nearzero || abs(walltime[j] - 2.0 * pi) < nearzero
        walltime[j] = Inf
      end
    end

    ## time till particle j hits or crosses wall
    movetime, j = findmin(walltime)
    if movetime == 0.0
      error("walking length zero")
    elseif movetime == Inf
      movetime = pi
    end

    totaltime += movetime
    if totaltime >= tune.traveltime
      movetime -= totaltime - tune.traveltime
      flag = true
    else
      tune.wallhits += 1
    end

    ## move the particle a time mt
    tune.velocity[:] = a * cos(movetime) - b * sin(movetime)
    tune.position[:] = a * sin(movetime) + b * cos(movetime)

    flag && break

    tune.position[j] = 0

    S1 = (S + ones(n)) / 2.0
    S1[j] = 0.0
    S2 = (S + ones(n)) / 2.0
    S2[j] = 1.0

    v2_new = tune.velocity[j]^2 +
             sign(tune.velocity[j]) * 2.0 * (logf(S2) - logf(S1))
    if v2_new > 0.0
      tune.velocity[j] = sqrt(v2_new) * sign(tune.velocity[j])
      S[j] *= -1.0
      tune.wallcrosses += 1
    else
      tune.velocity[j] *= -1.0
    end
  end

  ## convert from (-/+1) to (0/1)
  v[:] = (sign.(tune.position) + ones(n)) / 2.0
  v
end
