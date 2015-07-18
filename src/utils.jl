#################### Utility Macros ####################

macro modelexpr(args...)
  quote
    n = length($args)
    ex = Array(Expr, n)
    for i in 1:(n-1)
      x = $args[i]
      ex[i] = Expr(:(=), x, Expr(:ref, :model, QuoteNode(x)))
    end
    ex[n] = $args[n]
    Expr(:block, ex...)
  end
end

macro promote_scalarvariate(V)
  quote
    Base.promote_rule{T<:Real}(::Type{$V}, ::Type{T}) = Float64
  end
end


#################### Utility Functions ####################

dot(x) = dot(x, x)

invlogit(x) = 1.0 ./ (exp(-x) + 1.0)

logit(x) = log(x ./ (1.0 - x))

## Csorgo S and Faraway JJ. The exact and asymptotic distributions of the
## Cramer-von Mises statistic. Journal of the Royal Statistical Society,
## Series B, 58: 221-234, 1996.
function pcramer(q::Real)
  p = 0.0
  for k in 0:3
    c1 = 4.0 * k + 1.0
    c2 = c1^2 / (16.0 * q)
    p += gamma(k + 0.5) / factorial(k) * sqrt(c1) * exp(-c2) * besselk(0.25, c2)
  end
  p / (pi^1.5 * sqrt(q))
end
