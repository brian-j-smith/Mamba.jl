############ Binary Metropolized Gibbs Sampler ##############

#################### Sampler Constructor ####################

function BMG(params::Vector{Symbol})
  Sampler(params,
    quote
      x = round(Int,unlist(model,block,false)) #current state
      tunepar = tune(model, block)
      f = y -> logpdf!(model,y,block,false)
      bmg!(x,f)
      relist(model,x,block,false)
    end
  )
end

#################### Sampling Functions ####################

function bmg!(x::Vector{Int}, logf::Function)
  #proposal, initialize to current state
  y = x[:]
 
  m = fill(0.5,length(x))
  @inbounds for i in 1:length(x)
    x1 = x[:]; x0 = x[:]
    x1[i] = 1; x0[i] = 0

    x1l = logf(x1)
    x0l = logf(x0)

    # next line is equivalent to 
    # m[i] = exp(x1l)/(exp(x0l)+exp(x1l))
    # but more precise for log posteriors close to 1
    m[i] = 1/(2 + expm1(x0l-x1l))

    y[i] = rand(Bernoulli(m[i]))
  end
  if rand() < exp(logf(y) - logf(x))
    x[:] = y
  end
end
