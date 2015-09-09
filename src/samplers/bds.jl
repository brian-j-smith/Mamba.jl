##################### Binary Deterministic Sampler ##################

######################### Sampler Constructor #######################

function BDS{T <: Integer}(params::Vector{Symbol}, d::T, 
                           Γ::Vector{Vector{T}}=[[i] for i in 1:d])
  # If user supplies empty Γ, 
  # sample number of components for index set
  if size(Γ,1) == 0
    G = Truncated(Geometric(2/d), 1, d)
    k = rand(G)
    Γ = collect(combinations(1:d, k))
  end

  Sampler(params,
  quote
    x = round(Int,unlist(model,block,true))
    tunepar = tune(model,block)
    f = y -> logpdf!(model,y,block,true)
    bds!(x,tunepar["Γ"],f)
    relist(model,x,block,true)
  end,
  Dict("Γ" => Γ)
  )
end

########################## Sampling Functions #######################

function bds!{T <: Integer}(x::Vector{T}, Γ::Vector{Vector{T}},
                            logf::Function)
  #proposal, initialize to current state
  y = x[:]

  #indices to flip
  idx = Γ[rand(1:length(Γ))]

  y[idx] = 1 - x[idx]

  lr = logf(y) - logf(x)
  if rand() < exp(lr)
    x[:] = y
  end
end
