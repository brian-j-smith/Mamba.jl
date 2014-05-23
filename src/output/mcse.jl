function mcse{T<:Real}(x::Vector{T}, method::Symbol=:imse; args...)
  method == :imse ? mcse_imse(x) :
  method == :ipse ? mcse_ipse(x) :
  method == :bm ? mcse_bm(x; args...) :
    error("unsupported mcse method $method")
end

function mcse_bm{T<:Real}(x::Vector{T}; size::Integer=100)
  m = div(length(x), size)
  m >= 2 || error("2 or more batches needed to compute mcse")
  mbar = [mean(x[i * size + (1:size)]) for i in 0:m-1]
  sem(mbar)
end

function mcse_imse{T<:Real}(x::Vector{T})
  n = length(x)
  m = ifloor((n - 2.0) / 2.0)
  ghat = autocov(x, [0, 1])
  Ghat = sum(ghat)
  value = -ghat[1] + 2 * Ghat
  for i in 1:m
    Ghat = min(Ghat, sum(autocov(x, [2 * i, 2 * i + 1])))
    Ghat > 0 || break
    value += 2 * Ghat
  end
  sqrt(value / n)
end

function mcse_ipse{T<:Real}(x::Vector{T})
  n = length(x)
  m = ifloor((n - 2.0) / 2.0)
  ghat = autocov(x, [0, 1])
  value = ghat[1] + 2 * ghat[2]
  for i in 1:m
    Ghat = sum(autocov(x, [2 * i, 2 * i + 1]))
    Ghat > 0 || break
    value += 2 * Ghat
  end
  sqrt(value / n)
end
