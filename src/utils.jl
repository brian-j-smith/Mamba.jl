#################### Utility Functions ####################

centralrule(x) = cbrt(eps(eltype(x))) * max(one(eltype(x)), abs(x))

function gradient{T<:Real}(f::Function, x::Vector{T}, args...)
  n = length(x)
  g = Array(Float64, n)
  for i in 1:n
    epsilon = centralrule(x[i])
    x0 = x[i]
    x[i] = x0 - epsilon
    fx1 = f(x, args...)
    x[i] = x0 + epsilon
    fx2 = f(x, args...)
    x[i] = x0
    g[i] = (fx2 - fx1) / (epsilon + epsilon)
  end
  g
end

invlogit(x) = 1.0 / (exp(-x) + 1.0)

logit(x) = log(x ./ (1.0 - x))
