################################################################################
## Linear Regression
##   y ~ N(b0 + b1 * x, s2)
##   b0, b1 ~ N(0, 1000)
##   s2 ~ invgamma(0.001, 0.001)
################################################################################

using Mamba

## Data
data = [
3.7 0.825 1
3.5 1.09 1
1.25 2.5 1
0.75 1.5 1
0.8 3.2 1
0.7 3.5 1
0.6 0.75 0
1.1 1.7 0
0.9 0.75 0
0.9 0.45 0
0.8 0.57 0
0.55 2.75 0
0.6 3 0
1.4 2.33 1
0.75 3.75 1
2.3 1.64 1
3.2 1.6 1
0.85 1.415 1
1.7 1.06 0
1.8 1.8 1
0.4 2 0
0.95 1.36 0
1.35 1.35 0
1.5 1.36 0
1.6 1.78 1
0.6 1.5 0
1.8 1.5 1
0.95 1.9 0
1.9 0.95 1
1.6 0.4 0
2.7 0.75 1
2.35 0.03 0
1.1 1.83 0
1.1 2.2 1
1.2 2 1
0.8 3.33 1
0.95 1.9 0
0.75 1.9 0
1.3 1.625 1
]

vaso = Dict{Symbol, Any}(
  :y => data[:, end],
  :X => [ones(size(data,1)) mapslices(x -> (x - mean(x)) / sqrt(var(x)), 
                                      data[:, 1:(end-1)], 1)],
  :s2 => 10
)

## Log-transformed posterior + constant, gradient, hessian, and hessian 
logfgradhess = function(beta::DenseVector)
  eta = vaso[:X]*beta
  log_phi_eta = logcdf(Normal(),eta)
  log_phi_neg_eta = logcdf(Normal(), -eta)

  log_f_eta = logpdf(Normal(), eta)
  v=(vaso[:X]*beta.*(vaso[:y].*exp(log_f_eta - log_phi_eta)))

  grad = vaso[:X]' * (vaso[:y] .* exp(log_f_eta - log_phi_eta)) - 
         vaso[:X]' * ((1-vaso[:y]) .* exp(log_f_eta - log_phi_neg_eta)) - 
         beta / vaso[:s2]
  hess = vaso[:X]' * repmat(v', length(beta), 1) * vaso[:X] + 
         eye(length(beta)) / vaso[:s2]
  
  (log_f_eta, grad, hess)
end

## MCMC Simulation with Simplified Manifold 
## Metropolis-Adjusted Langevin Algorithm
## Without (1) and with (2) a user-specified proposal covariance matrix
n = 5000
sim1 = Chains(n, 2, names = ["b0", "b1"])
sim2 = Chains(n, 2, names = ["b0", "b1"])
epsilon = 0.1
Sigma = eye(2)
theta1 = SMMALAVariate([0.0, 0.0], epsilon, logfgradhess)
theta2 = SMMALAVariate([0.0, 0.0], epsilon, Sigma, logfgradhess)
for i in 1:n
  sample!(theta1)
  sample!(theta2)
  sim1[i, :, 1] = theta1
  sim2[i, :, 1] = theta2
end
describe(sim1)
describe(sim2)
