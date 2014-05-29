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


#################### Utility Functions ####################

dot(x) = dot(x, x)

invlogit(x) = 1.0 ./ (exp(-x) .+ 1.0)

logit(x) = log(x ./ (1.0 .- x))
