function gewekediag{T<:Real}(x::Vector{T}; first::Real=0.1, last::Real=0.5,
           etype=:bm, args...)
  n = length(x)
  x1 = x[1:int(first * n)]
  x2 = x[int(n - last * n + 1):n]
  z = (mean(x1) - mean(x2)) /
      sqrt(mcse(x, etype; args...)^2 + mcse(x, etype; args...)^2)
  [z, 1.0 - erf(abs(z) / sqrt(2.0))]
end

function gewekediag(c::MCMCChains; first::Real=0.1, last::Real=0.5,
           etype=:bm, args...)
  _, p, m = size(c)
  vals = Array(Float64, p, 2, m)
  for j in 1:p, k in 1:m
    vals[j,:,k] = gewekediag(c.value[:,j,k], first=first, last=last;
                    etype=etype, args...)
  end
  ChainSummary(vals, c.names, ["Z-score", "p-value"], header(c))
end
