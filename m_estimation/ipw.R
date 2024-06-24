#source('m_estimation/bias_correction.R')

ee.ipw <- function(data, known_e) {
  # Uses known propensity score which is assumed to be constant
  
  S <- data$S
  A <- data$A
  Y <- data$Y
  
  function(theta) {
    q  <- theta[1]
    ate_trial <- (theta[2] - theta[3])
    c(S - q,
      S / q * (A / known_e * Y - theta[2]),
      S / q * ((1 - A) / (1 - known_e) * Y - theta[3]),
      ate_trial - theta[4])
  }
}

estimate_ipw <- function(data, known_e = 1 / 2) {
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
    outer_args = list(known_e = known_e),
    corrections = list(bias_correction_.fay = correction(fay_bias_correction))
  )
  
  estimate = roots(out)[nparms]
  corrected_se = sqrt(get_corrections(out)$bias_correction_.fay[nparms, nparms])
  
  return(
    list(
      est =  estimate,
      se = corrected_se,
      ci.ll = estimate - 1.96 * corrected_se,
      ci.ul = estimate + 1.96 * corrected_se
    )
  )
  
  
  
}
