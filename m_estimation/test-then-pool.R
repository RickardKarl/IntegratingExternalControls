source('m_estimation/aipw.R')
source('m_estimation/pooling.R')
source('tests/lr_test.R')



test_then_pool <- function(data,
                           outcome_formula,
                           treatment_formula,
                           participant_formula,
                           propensity_score,
                           alpha = 0.05) {
  pval <- likelihood_ratio_test(data)
  
  
  if (pval < alpha) {
    out <-
      estimate_aipw(data,
                          outcome_formula,
                          treatment_formula,
                          propensity_score)
    
  } else {
    out <-
      estimate_pooling(
        data,
        outcome_formula,
        treatment_formula,
        participant_formula,
        propensity_score
      )
  }
  
  return(data.frame(out, reject_null = (pval < alpha)))
  
}