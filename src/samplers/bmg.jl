############ Binary Metropolized Gibbs Sampler ##############

#################### Types ####################

type BMGVariate <: VectorVariate
  value::Vector{Float64}

  BMGVariate(x::Vector{Float64}) = new(x)
end

#################### Sampler Constructor ####################

function BMG(params::Vector{Symbol})
  Sampler(params,
    quote
      x = unlist(model,block)
      f = y -> logpdf!(model,y,block)
      v = BMGVariate(x)
      bmg!(v,f)
      relist(model,v.value,block,false)
    end
  )
end

#################### Sampling Functions ####################

function bmg!(v::BMGVariate, logf::Function)
  #proposal, initialize to current state
  y = v[:]
 
  m = fill(0.5,length(v.value))
  for i in 1:length(v.value)
    v1 = v[:]; v0 = v[:]
    v1[i] = 1; v0[i] = 0

    v1l = logf(v1)
    v0l = logf(v0)

    m[i] = exp(v1l)/(exp(v0l)+exp(v1l))
    y[i] = rand(Bernoulli(m[i]))
  end
  if rand() < exp(logf(y) - logf(v.value))
    v[:] = y
  end
end
