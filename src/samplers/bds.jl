#################### Binary Deterministic Sampler ####################

#################### Types ####################

type BDSTune
  Γ::Vector{Vector{Int}}
end

type BDSVariate <: VectorVariate
  value::Vector{Float64}
  tune::BDSTune

  BDSVariate(x::Vector{Float64}, tune::BDSTune) = new(x, tune)
end

function BDSVariate(x::Vector{Float64}, tune=nothing)
  tune = BDSTune(
    Vector{Vector{Int}}[]
  )
  BDSVariate(x, tune)
end


#################### Sampler Constructor ####################

function BDS(params::Vector{Symbol}, d::Integer, k::Integer=1)
  d >= k > 0 || throw(ArgumentError("values must be d >= k > 0"))
  Γ = collect(combinations(1:d, k))
  BDS(params, Γ)
end

function BDS(params::Vector{Symbol}, Γ::Vector{Vector{Int}})
  Sampler(params,
    quote
      x = unlist(model, block, false)
      tunepar = tune(model, block)
      f = y -> logpdf!(model, y, block, false)
      v = BDSVariate(x)
      bds!(v, tunepar["Γ"], f)
      relist(model, v.value, block, false)
    end,
    Dict("Γ" => Γ)
  )
end


#################### Sampling Functions ####################

function bds!(v::BDSVariate, Γ::Vector{Vector{Int}}, logf::Function)
  x = v[:]

  # sample indices and flip values
  idx = Γ[rand(1:length(Γ))]
  x[idx] = 1.0 - v[idx]

  if rand() < exp(logf(x) - logf(v.value))
    v[:] = x
  end
  v.tune.Γ = Γ

  v
end
