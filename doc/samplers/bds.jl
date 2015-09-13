using Mamba

p = 10
n = 150

X = reshape(rand(Normal(0,1),p*n),n,p)
beta0 = rand(Normal(0,1),p)
gamma0 = rand(Bernoulli(0.5),p)

y = X*(beta0.*gamma0)

data = Dict(
  :y => y,
  :X => X,
  :beta0 => beta0,
  :n => n,
  :p => p
)

logf = function(x)
  logpdf(MvNormal(data[:X]*(data[:beta0].*x),1),data[:y])
end

t = 10000
gamma_names = fill("",p)
for i in 1:p
  gamma_names[i] = "gamma$i"
end
sim = Chains(t, p, names=gamma_names)
nu = BDSVariate(zeros(p))

for i in 1:t 
  bds!(nu,Array{Int64,1}[[j] for j in 1:p],logf)
  sim[i,:,1] = nu[1:p]
end

p=plot(sim,[:trace,:mixeddensity])
draw(p,filename="example")
describe(sim)