using Mamba
using Test
using Random

Random.seed!(123)

function simulate_NDARMA(N::Int, p::Int, q::Int, prob::Vector{Float64}, 
                         phi::Vector{Float64})
  X = zeros(Int64, N)
  X[1:p] = rand(Categorical(prob), p)
  d1 = Multinomial(1, phi)
  d2 = Categorical(prob)
  for t in (p+1):N
    alphabeta = rand(d1)
    eps = rand(d2, q + 1)
    X[t] = sum([X[(t-p):(t-1)]; eps] .* alphabeta)
  end
  return X
end

n = 10000
x1 = simulate_NDARMA(n, 1, 0, [0.25, 0.5, 0.05, 0.2], [0.8, 0.2])
x2 = simulate_NDARMA(n, 1, 0, [0.25, 0.5, 0.05, 0.2], [0.8, 0.2])
x3 = simulate_NDARMA(n, 1, 0, [0.25, 0.5, 0.05, 0.2], [0.8, 0.2])

sim = Chains(cat(x1,x2,x3,dims=3))
discretediag(sim)

res = Mamba.weiss(hcat(x1,x2,x3))
@test res[3] > 0.05


x1 = simulate_NDARMA(n, 1, 0, [0.25, 0.5, 0.05, 0.2], [0.8, 0.2])
x2 = simulate_NDARMA(n, 1, 0, [0.75, 0.15, 0.05, 0.05], [0.8, 0.2])
x3 = simulate_NDARMA(n, 1, 0, [0.25, 0.5, 0.05, 0.2], [0.8, 0.2])
res = Mamba.weiss(hcat(x1,x2,x3))
@test res[3] < 0.05
