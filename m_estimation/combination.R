combine_estimators <- function(est1, est2, var1, var2, cov) {
  lambda <- (var1 - cov) / (var1 + var2 - 2 * cov)
  est_combined <- (1 - lambda) * est1 + lambda * est2
  var_combined <- (var1 * var2 - cov ^ 2) / (var1 + var2 - 2 * cov)
  return(list(
    est = est_combined,
    var = var_combined,
    lambda = lambda
  ))
}
