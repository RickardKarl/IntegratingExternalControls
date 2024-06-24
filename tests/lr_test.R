library(lmtest)

likelihood_ratio_test <- function(data, outcome_formula, return_test_statistic = FALSE) {
  
  # Function to test if E[Y| X, S] = E[Y| X] using LR test
  
  data <- data[data$A==0,]
  
  full_model <- lm(Y ~ (.), data=data)
  
  subset <- dplyr::select(data,-c('S'))
  reduced_model <- lm(Y ~ (.), data=subset)
  
  out <- lrtest(full_model, reduced_model)
  if (!return_test_statistic){
    return(out$`Pr(>Chisq)`[2])
  } else{
    return(out$Chisq[2])
  }
}
