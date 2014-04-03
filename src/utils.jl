#################### Utility Functions ####################

invlogit(x) = 1.0 / (exp(-x) + 1.0)

logit(x) = log(x ./ (1.0 - x))
