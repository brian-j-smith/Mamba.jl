#################### Utility Functions ####################

function annotate(x::Matrix, rownames::Vector, colnames::Vector)
  hcat(["", rownames], vcat(colnames', x))
end

invlogit(x) = 1.0 / (exp(-x) + 1.0)

logit(x) = log(x ./ (1.0 - x))
