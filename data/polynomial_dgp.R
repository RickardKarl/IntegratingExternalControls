# A data-generating process
library(MASS)

baseline <- function(x) {
  return(
    1 / 2 * x$X1 + x$X2 - 1 / 2 * x$X3 +  x$X4 - 1 / 2 * x$X5
    - 1 / 4 * x$X1 ^ 2 - x$X2 ^ 2 - 1 / 2 * x$X3 ^ 2 -  x$X4 ^ 2 - 1 / 2 * x$X5 ^ 2
    + 1 / 2 * x$X6 ^ 2 + 1 / 2 * x$X7 ^ 2 + 1 / 2 * x$X8 ^ 2  + 1 / 2 * x$X9 ^ 2 +  1 / 2 * x$X10 ^ 2
  )
}

modifier <- function(x) {
  return(5)
}

dgp <- function(n1, n0, delta) {
  # n1 : samples from trial, S=1
  # n0 : samples from target population, S=0
  # beta: controls violation of E[Y^0 \mid X, S=1] = E[Y^0 \mid X, S=0]
  
  n = n1 + n0
  # Study assignment
  S0 <- rep.int(0, n0)
  S1 <- rep.int(1, n1)
  S <- c(S0, S1)
  
  # Treatment assignment
  A0_S0 <- rep.int(0, length(S0))
  A0_S1 <- rep.int(0, length(S1) %/% 2)
  A1_S1 <- rep.int(1, n1 - length(A0_S1))
  A <- c(A0_S0, A0_S1, A1_S1)
  
  # Covariates
  sigma = matrix(rep(0, (10) * (10)), nrow = (10))
  diag(sigma) = 1
  X_S0 <- mvrnorm(n = n0, rep.int(delta, 10),  sigma)
  X_S1 <- mvrnorm(n = n1, rep.int(0, 10),  sigma)
  X = rbind(data.frame(X_S0), data.frame(X_S1))
  
  # Outcome
  Y <- baseline(X) + A * modifier(X) + rnorm(n, sd = 1)
  
  return(data.frame(A, Y, X, S))
}

true_ate <- function() {
  return(5)
}

true_control_mean <- function() {
  return(-.75) # approximate
}
