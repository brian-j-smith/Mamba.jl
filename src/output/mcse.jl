function mcse{T<:Real}(x::Vector{T}, method::Symbol=:bm; args...)
  method == :bm ? mcse_bm(x; args...) :
    error("unsupported mcse method $method")
end

function mcse_bm{T<:Real}(x::Vector{T}; size::Integer=100)
  m = div(length(x), size)
  m >= 2 || error("2 or more batches needed to compute mcse")
  mbar = [mean(x[i * size + (1:size)]) for i in 0:m-1]
  sem(mbar)
end
