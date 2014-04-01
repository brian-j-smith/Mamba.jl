#################### Gelman, Rubin, and Brooks Diagnostics ####################

function gelmandiag(c::MCMCChain; alpha::Real=0.05, mpsrf::Bool=false,
                    transform::Bool=false)
  n, p, m = size(c.data)
  m >= 2 || error("2 or more chains needed to run gelman diagnostic")

  psi = transform ? transform_support(c) : c.data

  S2 = mapslices(x -> reshape(crosscov(x, x, [0]), p, p), psi, [1,2])
  W = mapreduce(i -> S2[:,:,i], +, 1:m) / m

  psibar = reshape(mapslices(mean, psi, 1), p, m)'
  B = n * reshape(crosscov(psibar, psibar, [0]), p, p)

  w = diag(W)
  b = diag(B)
  s2 = reshape(mapslices(diag, S2, [1,2]), p, m)'
  psibar2 = map(i -> mean(psibar[:,i]), 1:p)

  var_w = map(i -> var(s2[:,i]), 1:p) / m
  var_b = (2 / (m - 1)) * b.^2
  var_wb = (n / m) * (diag(reshape(crosscov(s2, psibar.^2, [0]), p, p)) -
           2 * psibar2 .* diag(reshape(crosscov(s2, psibar, [0]), p, p)))

  V = ((n - 1) / n) * w + ((m + 1) / (m * n)) * b
  var_V = ((n - 1)^2 * var_w + ((m + 1) / m)^2 * var_b +
           (2 * (n - 1) * (m + 1) / m) * var_wb) / n^2
  df = 2 * V.^2 ./ var_V
  B_df = m - 1
  W_df = 2 * w.^2 ./ var_w

  R_fixed = (n - 1) / n
  R_random = ((m + 1) / (m * n)) * b ./ w
  R_est = R_fixed + R_random
  q = 1 - alpha / 2
  R_upper = R_fixed + R_random .*
            map(df2 -> quantile(FDist(B_df, df2), q), W_df)

  psrf = sqrt((df + 3) / (df + 1) * [R_est R_upper])
  psrf_labels = ["PSRF", string(100 * q) * "%"]
  psrf_names = c.names

  if mpsrf
    x = [(n - 1) / n + (m + 1) / m * eigmax(inv(cholfact(W)) * B) / n NaN]
    psrf = vcat(psrf, x)
    psrf_names = [psrf_names, "Multivariate"]
  end

  ChainSummary(psrf, psrf_names, psrf_labels, header(c))
end
