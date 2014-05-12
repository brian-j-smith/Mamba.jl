using MCMCsim
using Distributions

## Data
leuk = (String => Any)[
  "t_obs" => [1, 1, 2, 2, 3, 4, 4, 5, 5, 8, 8, 8, 8, 11, 11, 12, 12, 15, 17, 22,
              23, 6, 6, 6, 6, 7, 9, 10, 10, 11, 13, 16, 17, 19, 20, 22, 23, 25,
              32, 32, 34, 35],
  "fail" => [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
             1, 1, 0, 1, 0, 1, 0, 0, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0],
  "Z" => [0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5,
          0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5,
          -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5,
          -0.5, -0.5, -0.5, -0.5],
  "t" => [1, 2, 3, 4, 5, 6, 7, 8, 10, 11, 12, 13, 15, 16, 17, 22, 23, 35]
]
leuk["N"] = N = length(leuk["t_obs"])
leuk["T"] = T = length(leuk["t"]) - 1

leuk["Y"] = Array(Integer, N, T)
leuk["dN"] = Array(Integer, N, T)
for i in 1:N
  for j in 1:T
    leuk["dN"][i,j] = leuk["fail"][i] * int(leuk["t_obs"][i] == leuk["t"][j])
    leuk["Y"][i,j] = int(leuk["t_obs"][i] >= leuk["t"][j])
  end
end

leuk["c"] = 0.001
leuk["r"] = 0.1


## Model Specification

model = MCMCModel(

  dN = MCMCStochastic(2,
    @modelexpr(Y, beta, Z, dL0, N, T,
      Distribution[
        Y[i,j] > 0 ? Poisson(exp(beta * Z[i]) * dL0[j]) : Flat()
        for i in 1:N, j in 1:T
      ]
    ),
    false
  ),

  mu = MCMCLogical(1,
    @modelexpr(c, r, t,
      c * r * (t[2:end] - t[1:end-1])
    ),
    false
  ),

  dL0 = MCMCStochastic(1,
    @modelexpr(mu, c, T,
      Distribution[Gamma(mu[j], c) for j in 1:T]
    ),
    false
  ),

  beta = MCMCStochastic(
    :(Normal(0, 1000))
  )

)


## Initial Values
inits = [
  ["dN" => leuk["dN"], "beta" => 0, "dL0" => ones(leuk["T"])],
  ["dN" => leuk["dN"], "beta" => 1, "dL0" => fill(2, leuk["T"])]
]


## Sampling Scheme
scheme = [AMWG(["dL0"], fill(0.01, leuk["T"])),
          Slice(["beta"], [3.0])]
setsamplers!(model, scheme)


## MCMC Simulations
sim = mcmc(model, leuk, inits, 10000, burnin=2500, thin=2, chains=2)
describe(sim)
