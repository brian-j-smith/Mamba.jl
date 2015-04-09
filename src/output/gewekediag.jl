function gewekediag{T<:Real}(x::Vector{T}; first::Real=0.1, last::Real=0.5,
           etype=:imse, args...)
  if !(0 < first < 1)
    error("first must be in (0, 1)")
  elseif !(0 < last < 1)
    error("last must be in (0, 1)")
  elseif first + last > 1
    error("first and last sequences are overlapping")
  end
  n = length(x)
  x1 = x[1:int(first * n)]
  x2 = x[int(n - last * n + 1):n]
  z = (mean(x1) - mean(x2)) /
      sqrt(mcse(x1, etype; args...)^2 + mcse(x2, etype; args...)^2)
  [round(z, 3), round(1.0 - erf(abs(z) / sqrt(2.0)), 4)]
end

function gewekediag(c::Chains; first::Real=0.1, last::Real=0.5, etype=:imse,
           args...)
  _, p, m = size(c.value)
  vals = Array(Float64, p, 2, m)
  for j in 1:p, k in 1:m
    vals[j,:,k] = gewekediag(c.value[:,j,k], first=first, last=last;
                    etype=etype, args...)
  end
  ChainSummary(vals, c.names, ["Z-score", "p-value"], header(c))
end
