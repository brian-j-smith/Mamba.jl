using MCMCsim
using Distributions

## Data
inhalers = (String => Any)[
  "pattern" =>
    [1 1 1 1 2 2 2 2 3 3 3 3 4 4 4 4
     1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4]',
  "Ncum" =>
    [ 59 157 173 175 186 253 270 271 271 278 280 281 282 285 285 286
     122 170 173 175 226 268 270 271 278 280 281 281 284 285 286 286]',
  "treat" =>
    [ 1 -1
     -1  1],
  "period" =>
    [1 -1
     1 -1],
  "carry" =>
    [0 -1
     0  1],
  "N" => 286,
  "T" => 2,
  "G" => 2,
  "Npattern" => 16,
  "Ncut" => 3
]

inhalers["group"] = group = Array(Integer, inhalers["N"])
inhalers["response"] = response = Array(Integer, inhalers["N"], inhalers["T"])

for i in 1:inhalers["Ncum"][1,1]
  group[i] = 1
  for t in 1:inhalers["T"]
    response[i,t] = inhalers["pattern"][1,t]
  end
end

for i in inhalers["Ncum"][1,1]+1 : inhalers["Ncum"][1,2]
  group[i] = 2
  for t in 1:inhalers["T"]
    response[i,t] = inhalers["pattern"][1,t]
  end
end

for k in 2:inhalers["Npattern"]
  for i in inhalers["Ncum"][k-1,2]+1 : inhalers["Ncum"][k,1]
    group[i] = 1
    for t in 1:inhalers["T"]
      response[i,t] = inhalers["pattern"][k,t]
    end
  end
  for i in inhalers["Ncum"][k,1]+1 : inhalers["Ncum"][k,2]
    group[i] = 2
    for t in 1:inhalers["T"]
      response[i,t] = inhalers["pattern"][k,t]
    end
  end
end


## Model Specification

model = MCMCModel(

  response = MCMCStochastic(2,
    @modelexpr(a1, a2, a3, mu, group, b, N, T,
      begin
        a = Float64[a1, a2, a3]
        p = Array(Float64, 4)
        Distribution[
          begin
            p[1] = 1.0
            for j in 1:3
              Q = invlogit(-(a[j] + mu[group[i],t] + b[i]))
              p[j] -= Q
              p[j+1] = Q
            end
            Categorical(p)
          end
          for i in 1:N, t in 1:T
        ]
      end
    ),
    false
  ),

  mu = MCMCLogical(2,
    @modelexpr(beta, treat, pi, period, kappa, carry, G, T,
      [
        beta * treat[g,t] / 2 + pi * period[g,t] / 2 + kappa * carry[g,t]
        for g in 1:G, t in 1:T
      ]
    ),
    false
  ),

  b = MCMCStochastic(1,
    @modelexpr(s2, N,
      IsoNormal(N, sqrt(s2))
    ),
    false
  ),

  a1 = MCMCStochastic(
    @modelexpr(a2,
      Uniform(-1000, a2)
    )
  ),

  a2 = MCMCStochastic(
    @modelexpr(a3,
      Uniform(-1000, a3)
    )
  ),

  a3 = MCMCStochastic(
    :(Uniform(-1000, 1000))
  ),

  beta = MCMCStochastic(
    :(Normal(0, 1000))
  ),

  pi = MCMCStochastic(
    :(Normal(0, 1000))
  ),

  kappa = MCMCStochastic(
    :(Normal(0, 1000))
  ),

  s2 = MCMCStochastic(
    :(InverseGamma(0.001, 0.001))
  )

)


## Initial Values
inits = [
  ["response" => inhalers["response"], "beta" => 0, "pi" => 0, "kappa" => 0,
   "a1" => 2, "a2" => 3, "a3" => 4, "s2" => 1, "b" => ones(inhalers["N"])],
  ["response" => inhalers["response"], "beta" => 1, "pi" => 1, "kappa" => 0,
   "a1" => 3, "a2" => 4, "a3" => 5, "s2" => 10, "b" => ones(inhalers["N"])]
]


## Sampling Scheme
scheme = [AMWG(["b"], fill(0.1, inhalers["N"])),
          Slice(["a1", "a2", "a3"], fill(1.0, 3)),
          AMM(["beta", "pi", "kappa"], 0.1 * eye(3)),
          Slice(["s2"], [1.0])]
setsamplers!(model, scheme)


## MCMC Simulations
sim = mcmc(model, inhalers, inits, 5000, burnin=1000, thin=2, chains=2)
describe(sim)
