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
  sigma = matrix(rep(1 / 10, (10) * (10)), nrow = (10))
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

# 
# data <- dgp(100, 200, 1/2)
# hist(data$X1[data$S == 1], breaks = 20)
# hist(data$X1[data$S == 0], breaks = 20)
# # 
# # Fit misspecified regression model
#  
# print(paste('correct', mean(data$Y[(data$S == 1) &
#                                       (data$A == 1)]) - mean(data$Y[(data$S == 1) &
#                                                                       (data$A == 0)])))
# print(paste('incorrect', mean(data$Y[(data$S == 1) &
#                                         (data$A == 1)]) - mean(data$Y[(data$A == 0)])))
 
#S0.model <- lm(Y~X1+X2+X3+X4+X5+X6+X7+X8+X9+X10, data, subset=())

# 
# data <- dgp(1000, 1000, 1/2)
# hist(data$X1[data$S == 1], breaks = 20)
# hist(data$X1[data$S == 0], breaks = 20)
# # Fit linear regression for S~X
# partici_model <-
#   glm(S ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10,
#       data,
#       family = binomial)
# summary(partici_model)
# weights <- predict(partici_model, newdata = data, type = "response")
# hist(weights, breaks=40)

#
#true_ate()
#
#
# Check that IPW is stable
# source('estimators/ipw.R')
# iterations <- 1000
# for (j in 1:10) {
#   estimates <- matrix(nrow = iterations)
#   for (i in 1:iterations) {
#     data <- dgp(100, 100, 0)
#     estimates[i] <- ipw.ate(data)$est
#   }
#   cat(j, mean(estimates) - 5, var(estimates), '\n')
# }

#
#
#
# Check that h model has better fit
# g_error <- matrix(nrow = 25)
# h_error <- matrix(nrow = 25)
# for (i in 1:25) {
#   data <- dgp(100, 200, 1/2)
#   outcome_formula <-
#     formula(
#       Y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10
#     )
# 
#   # g_0 model
#   trial_model <- lm(outcome_formula, data, subset = (A == 0 & S==1))
#   # h model
#   partici_model <-
#     glm(S ~ X1 + X2 + X3 + X4 + X5,
#         data,
#         subset = (A == 0),
#         family = binomial)
#   ec_model <-
#     lm(
#       outcome_formula,
#       data,
#       subset = (A == 0),
#       weights = predict(partici_model, newdata = data, type = "response")
#     )
# 
#   # Compute residuals
#   test_data <- data <- dgp(5000, 500, 0)
#   trial_error <-
#     (predict(trial_model, newdata = test_data) - test_data$Y)[(test_data$A == 0) &
#                                                                 (test_data$S == 1)]
#   ec_error <-
#     (predict(ec_model, newdata = test_data) - test_data$Y)[(test_data$A == 0) &
#                                                              (test_data$S == 1)]
# 
#   g_error[i] <- mean((trial_error) ^ 2)
#   h_error[i] <- mean((ec_error) ^ 2)
#   #cat('g_0 model', mean((trial_error) ^ 2), '\n')
#   #cat('h model', mean((ec_error) ^ 2), '\n')
# }
# cat('g error', mean(g_error), var(g_error), '\n')
# cat('h error', mean(h_error), var(h_error), '\n')

#
# X <-
#   dplyr::select(data, c('X1', 'X2', 'X3', 'X4', 'X5'))[data$S == 1, ]
#
# p1 <- ggplot(data, aes(
#   x = X1,
#   y = Y,
#   color = A,
#   group = A
# )) +
#   geom_point()
# p2 <- ggplot(data, aes(
#   x = X2,
#   y = Y,
#   color = A,
#   group = A
# )) +
#   geom_point()
# p3 <- ggplot(data, aes(
#   x = X3,
#   y = Y,
#   color = A,
#   group = A
# )) +
#   geom_point()
#
# p4 <- ggplot(data, aes(
#   x = X4,
#   y = Y,
#   color = A,
#   group = A
# )) +
#   geom_point()
#
# p5 <- ggplot(data, aes(
#   x = X5,
#   y = Y,
#   color = A,
#   group = A
# )) +
#   geom_point()
# grid.arrange(p1, p2, p3, p4, p5)
#
# p1 <- ggplot(data, aes(
#   x = X6,
#   y = Y,
#   color = A,
#   group = A
# )) +
#   geom_point()
# p2 <- ggplot(data, aes(
#   x = X7,
#   y = Y,
#   color = A,
#   group = A
# )) +
#   geom_point()
# p3 <- ggplot(data, aes(
#   x = X8,
#   y = Y,
#   color = A,
#   group = A
# )) +
#   geom_point()
#
# p4 <- ggplot(data, aes(
#   x = X9,
#   y = Y,
#   color = A,
#   group = A
# )) +
#   geom_point()
#
# p5 <- ggplot(data, aes(
#   x = X10,
#   y = Y,
#   color = A,
#   group = A
# )) +
#   geom_point()
# grid.arrange(p1, p2, p3, p4, p5)
#
#