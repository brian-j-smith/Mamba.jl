using MCMCsim
using Distributions

## Data
data = (String => Any)[
  "y" => [151, 199, 246, 283, 320,
          145, 199, 249, 293, 354,
          147, 214, 263, 312, 328,
          155, 200, 237, 272, 297,
          135, 188, 230, 280, 323,
          159, 210, 252, 298, 331,
          141, 189, 231, 275, 305,
          159, 201, 248, 297, 338,
          177, 236, 285, 350, 376,
          134, 182, 220, 260, 296,
          160, 208, 261, 313, 352,
          143, 188, 220, 273, 314,
          154, 200, 244, 289, 325,
          171, 221, 270, 326, 358,
          163, 216, 242, 281, 312,
          160, 207, 248, 288, 324,
          142, 187, 234, 280, 316,
          156, 203, 243, 283, 317,
          157, 212, 259, 307, 336,
          152, 203, 246, 286, 321,
          154, 205, 253, 298, 334,
          139, 190, 225, 267, 302,
          146, 191, 229, 272, 302,
          157, 211, 250, 285, 323,
          132, 185, 237, 286, 331,
          160, 207, 257, 303, 345,
          169, 216, 261, 295, 333,
          157, 205, 248, 289, 316,
          137, 180, 219, 258, 291,
          153, 200, 244, 286, 324]
]
data["i"] = Integer[div(i - 1, 5) + 1 for i in 1:150]
data["j"] = Integer[(i - 1) % 5 + 1 for i in 1:150]
data["x"] = [8, 15, 22, 29, 36][data["j"]]
data["xm"] = data["x"] - 22.0


## Model Specification

rats = MCMCModel(

  s2_c = MCMCStochastic(
    :(InverseGamma(0.001, 0.001))
  ),

  mu_alpha = MCMCStochastic(
    :(Normal(0.0, 1.0e6)),
    false
  ),

  s2_alpha = MCMCStochastic(
    :(InverseGamma(0.001, 0.001)),
    false
  ),

  mu_beta = MCMCStochastic(
    :(Normal(0.0, 1.0e6))
  ),

  s2_beta = MCMCStochastic(
    :(InverseGamma(0.001, 0.001)),
    false
  ),

  alpha = MCMCStochastic(30,
    quote
      mu = model["mu_alpha"] * ones(30)
      s2 = model["s2_alpha"]
      IsoNormal(mu, sqrt(s2))
    end,
    false
  ),

  beta = MCMCStochastic(30,
    quote
      mu = model["mu_beta"] * ones(30)
      s2 = model["s2_beta"]
      IsoNormal(mu, sqrt(s2))
    end,
    false
  ),

  alpha0 = MCMCLogical(
    :(model["mu_alpha"] - 22.0 * model["mu_beta"])
  ),

  y = MCMCStochastic(150,
    quote
      alpha = model["alpha"][model["i"]]
      beta = model["beta"][model["i"]]
      mu = alpha + beta .* model["xm"]
      IsoNormal(mu, sqrt(model["s2_c"]))
    end,
    false
  )

)


## Initial Values
inits = [
  ["y" => data["y"], "alpha" => fill(250, 30), "beta" => fill(6, 30),
   "mu_alpha" => 150, "mu_beta" => 10, "s2_c" => 1, "s2_alpha" => 1,
   "s2_beta" => 1],
  ["y" => data["y"], "alpha" => fill(20, 30), "beta" => fill(0.6, 30),
   "mu_alpha" => 15, "mu_beta" => 1, "s2_c" => 10, "s2_alpha" => 10,
   "s2_beta" => 10]
]


## Sampling Scheme
scheme = [SamplerSlice(["s2_c"], [10.0]),
          SamplerAMWG(["alpha"], 100 * ones(30), adapt=:all),
          SamplerSliceWG(["mu_alpha", "s2_alpha"], [100.0, 10.0]),
          SamplerAMWG(["beta"], ones(30), adapt=:all),
          SamplerSliceWG(["mu_beta", "s2_beta"], [1.0, 1.0])]
setsamplers!(rats, scheme)


## MCMC Simulations
sim = mcmc(rats, data, inits, 10000, burnin=2500, thin=2, chains=2)
describe(sim)
