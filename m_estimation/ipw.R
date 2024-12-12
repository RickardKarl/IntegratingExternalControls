ee.ipw <- function(data, propensity_score) {
  S <- data$S
  A <- data$A
  Y <- data$Y
  
  function(theta) {
    q  <- theta[1]
    ate_trial <- (theta[2] - theta[3])
    c(S - q,
      S / q * (A / propensity_score * Y - theta[2]),
      S / q * ((1 - A) / (1 - propensity_score) * Y - theta[3]),
      ate_trial - theta[4])
  }
}

estimate_ipw <-
  function(data,
           propensity_score) {
    
    ############################################################################
    # Do m-estimation 
    ############################################################################
    
    require(geex)
    
    start_val <-
      c(mean(data$S),
        mean(data$Y[(data$A == 1) & (data$S == 1)]),
        mean(data$Y[(data$A == 0) & (data$S == 1)]),
        mean(data$Y[(data$A == 1) &
                      (data$S == 1)]) -  mean(data$Y[(data$A == 0) &
                                                       (data$S == 1)]))
    nparms <- length(start_val)
    
    out <- m_estimate(
      estFUN = ee.ipw,
      data = data,
      root_control = setup_root_control(start = start_val),
      outer_args = list(propensity_score = propensity_score),
      corrections =  list(bias_correction_.fay = correction(fay_bias_correction))
    )
    
    ############################################################################
    
    # get estimates of ATE
    estimate = roots(out)[nparms]
    se = sqrt(vcov(out)[nparms, nparms])
    corr_se = sqrt(get_corrections(out)$bias_correction_.fay[nparms, nparms])
    
    # get estimates of control mean outcome
    estimate.A0 = roots(out)[nparms - 1]
    se.A0 = sqrt(vcov(out)[nparms - 1, nparms - 1])
    corr_se.A0 = sqrt(get_corrections(out)$bias_correction_.fay[nparms - 1, nparms - 1])
    
    return(
      list(
        est =  estimate,
        se = se,
        ci.ll = estimate - 1.96 * se,
        ci.ul = estimate + 1.96 * se,
        corr_se = corr_se,
        corr_ci.ll = estimate - 1.96 * corr_se,
        corr_ci.ul = estimate + 1.96 * corr_se,
        est.A0 = estimate.A0,
        se.A0 = se.A0,
        ci.ll.A0 = estimate.A0 - 1.96 * se.A0,
        ci.ul.A0 = estimate.A0 + 1.96 * se.A0,
        corr_se.A0 = corr_se.A0,
        corr_ci.ll.A0 = estimate.A0 - 1.96 * corr_se.A0,
        corr_ci.ul.A0 = estimate.A0 + 1.96 * corr_se.A0
      )
    )
    
    
    
  }
