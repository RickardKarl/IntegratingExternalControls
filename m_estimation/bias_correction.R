
# Finite sample bias correction for variance estimation
bias_correction <- function(components, b) {
  A <- grab_bread(components)
  A_i <- grab_bread_list(components)
  B_i <- grab_meat_list(components)
  Ainv <- solve(A)
  
  H_i <- lapply(A_i, function(m) {
    diag((1 - pmin(b, diag(m %*% Ainv))) ^ (-0.5))
  })
  
  Bbc_i <- lapply(seq_along(B_i), function(i) {
    H_i[[i]] %*% B_i[[i]] %*% H_i[[i]]
  })
  Bbc   <- apply(simplify2array(Bbc_i), 1:2, sum)
  
  compute_sigma(A = A, B = Bbc)
  
}