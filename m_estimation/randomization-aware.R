#source('m_estimation/bias_correction.R')


ee.robin <- function(data, models) {
  # Get design matrix
  Xe  <-
    grab_design_matrix(data, rhs_formula = grab_fixed_formula(models$e))
  Xg0 <-
    grab_design_matrix(data, rhs_formula = grab_fixed_formula(models$g0))
  Xg1 <-
    grab_design_matrix(data, rhs_formula = grab_fixed_formula(models$g1))
  Xh  <-
    grab_design_matrix(data, rhs_formula = grab_fixed_formula(models$h))
  Xeta <-
    grab_design_matrix(data, rhs_formula = grab_fixed_formula(models$eta))
  
  # Get indices for parameters
  e_pos  <- 1:ncol(Xe)
  g0_pos <- (max(e_pos)  + 1):(max(e_pos)  + ncol(Xg0))
  g1_pos <- (max(g0_pos) + 1):(max(g0_pos) + ncol(Xg1))
  h_pos  <- (max(g1_pos) + 1):(max(g1_pos) + ncol(Xh))
  eta_pos <- (max(h_pos)  + 1):(max(h_pos)  + ncol(Xeta))
  
  # Grab estimating function from each model object
  e_scores  <- grab_psiFUN(models$e, data)
  g0_scores <- grab_psiFUN(models$g0, data)
  g1_scores <- grab_psiFUN(models$g1, data)
  h_scores <- grab_psiFUN(models$h, data)
  eta_scores <- grab_psiFUN(models$eta, data)
  
  S <- data$S
  A <- data$A
  Y <- data$Y
  
  function(theta) {
    e  <- plogis(Xe %*% theta[e_pos])[, 1]
    g0 <- Xg0 %*% theta[g0_pos]
    g1 <- Xg1 %*% theta[g1_pos]
    h  <- Xh %*% theta[h_pos]
    eta <- plogis(Xeta %*% theta[eta_pos])[, 1]
    q  <- theta[length(theta) - 5]
    ate_trial <-
      (theta[length(theta) - 4] - theta[length(theta) - 3])
    ate_robin <-
      (theta[length(theta) - 4] - theta[length(theta) - 2])
    c(
      S * e_scores(theta[e_pos]),
      S * (1 - A) * g0_scores(theta[g0_pos]),
      S * A * g1_scores(theta[g1_pos]),
      (1 - A) * eta * e / (1 - e) ^ 2 * h_scores(theta[h_pos]),
      (1 - A) * eta_scores(theta[eta_pos]),
      S - q,
      S / q * (A / e * (Y - g1) + g1 - theta[length(theta) - 4]),
      S / q * ((1 - A) / (1 - e) * (Y - g0) + g0 - theta[length(theta) - 3]),
      # no pooling
      S / q * ((1 - A) / (1 - e) * (Y - h) + h - theta[length(theta) - 2]),
      # robust pooling
      ate_trial - theta[length(theta) - 1],
      ate_robin - theta[length(theta)]
    )
  }
}


estimate_robin <-
  function(data,
           outcome_formula,
           treatment_formula,
           participant_formula) {
    require(geex)
    e_model  <-
      glm(
        treatment_formula,
        subset = (S == 1),
        data = data,
        family = binomial
      )
    eta_model <-
      glm(
        participant_formula,
        subset = (A == 0),
        data = data,
        family = binomial
      )
    
    g0_model <-
      glm(outcome_formula,
          subset = (A == 0 & S == 1),
          data = data)
    
    g1_model <-
      glm(outcome_formula,
          subset = (A == 1 & S == 1),
          data = data)
    
    h_model <-
      glm(outcome_formula,
          subset = (A == 0),
          data = data)
    
    
    models_robin <-
      list(
        e = e_model,
        g0 = g0_model,
        g1 = g1_model,
        h = h_model,
        eta = eta_model
      )
    nparms_robin <-
      sum(unlist(lapply(models_robin, function(x)
        length(coef(
          x
        ))))) + 6
    
    start_val <-
      c(
        coef(e_model),
        coef(g0_model),
        coef(g1_model),
        coef(h_model),
        coef(eta_model),
        mean(data$S),
        mean(data$Y[(data$A == 1) & (data$S == 1)]),
        mean(data$Y[(data$A == 0) & (data$S == 1)]),
        mean(data$Y[(data$A == 0) & (data$S == 1)]),
        mean(data$Y[(data$A == 1) &
                      (data$S == 1)]) -  mean(data$Y[(data$A == 0) &
                                                       (data$S == 1)]),
        mean(data$Y[(data$A == 1) &
                      (data$S == 1)]) -  mean(data$Y[(data$A == 0) &
                                                       (data$S == 1)])
      )
    
    out <- m_estimate(
      estFUN = ee.robin,
      data = data,
      root_control = setup_root_control(start = start_val),
      outer_args = list(models = models_robin),
      corrections = list(bias_correction_.fay = correction(fay_bias_correction))
    )
    
    # Get estimates
    est_trial <- roots(out)[nparms_robin - 1]
    est_robin <- roots(out)[nparms_robin]
    
    
    # Compute combined estimator
    uncorr_var_trial <-
      vcov(out)[nparms_robin - 1 , nparms_robin - 1]
    uncorr_var_robin <- vcov(out)[nparms_robin, nparms_robin]
    uncorr_cov_robin_trial <-
      vcov(out)[nparms_robin, nparms_robin - 1]
    
    lambda <-
      (uncorr_var_trial - uncorr_cov_robin_trial) / (uncorr_var_trial + uncorr_var_robin - 2 * uncorr_cov_robin_trial)
    est_combined <- (1 - lambda) * est_trial + lambda * est_robin
    
    
    
    # Compute bias-corrected SEs
    corr_var_trial <-
      get_corrections(out)$bias_correction_.fay[nparms_robin - 1, nparms_robin - 1]
    corr_var_robin <-
      get_corrections(out)$bias_correction_.fay[nparms_robin, nparms_robin]
    corr_cov_robin_trial <-
      get_corrections(out)$bias_correction_.fay[nparms_robin, nparms_robin - 1]
    corr_var_combined <-
      (corr_var_trial * corr_var_robin - corr_cov_robin_trial ^ 2) / (corr_var_trial + corr_var_robin - 2 * corr_cov_robin_trial)
    
    # Confidence intervals
    ci.ll_trial <- est_trial - 1.96 * sqrt(corr_var_trial)
    ci.ul_trial <- est_trial + 1.96 * sqrt(corr_var_trial)
    ci.ll_robin <- est_robin - 1.96 * sqrt(corr_var_robin)
    ci.ul_robin <- est_robin + 1.96 * sqrt(corr_var_robin)
    ci.ll_combined <- est_combined - 1.96 * sqrt(corr_var_combined)
    ci.ul_combined <- est_combined + 1.96 * sqrt(corr_var_combined)
    
    
    return(
      list(
        est_trial = est_trial,
        est_robin = est_robin,
        est_combined = est_combined,
        se_trial = sqrt(corr_var_trial),
        se_robin = sqrt(corr_var_robin),
        se_combined = sqrt(corr_var_combined),
        ci.ll_trial = ci.ll_trial,
        ci.ul_trial = ci.ul_trial,
        ci.ll_robin = ci.ll_robin,
        ci.ul_robin = ci.ul_robin,
        ci.ll_combined = ci.ll_combined,
        ci.ul_combined = ci.ul_combined
      )
    )
    
  }
