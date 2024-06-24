source('m_estimation/no_pooling.R')
source('m_estimation/pooling.R')
source('tests/lr_test.R')



test_then_pool <- function(data,
                           outcome_formula,
                           treatment_formula,
                           participant_formula,
                           alpha = 0.05) {
  pval <- likelihood_ratio_test(data)
  
  
  if (pval < 0.05) {
    out <- estimate_no_pooling(data, outcome_formula, treatment_formula)
    
  } else {
    out <-
      estimate_pooling(data,
                       outcome_formula,
                       treatment_formula,
                       participant_formula)
  }
  
  return(data.frame(out, reject_null = (pval < alpha)))
  
}