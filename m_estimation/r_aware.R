ee.r_aware <-
  function(data,
           models,
           propensity_score = NULL,
           use_no_weights_h = FALSE) {
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
      if (!is.null(propensity_score)) {
        e[] <- propensity_score
      }
      g0 <- Xg0 %*% theta[g0_pos]
      g1 <- Xg1 %*% theta[g1_pos]
      h  <- Xh %*% theta[h_pos]
      eta <- plogis(Xeta %*% theta[eta_pos])[, 1]
      h_score_weights <- eta * e / (1 - e) ^ 2
      if (use_no_weights_h) {
        h_score_weights[] <- 1.0
      }
      q  <- theta[length(theta) - 5]
      ate_trial <-
        (theta[length(theta) - 4] - theta[length(theta) - 3])
      ate_r_aware <-
        (theta[length(theta) - 4] - theta[length(theta) - 2])
      
      c(
        S * e_scores(theta[e_pos]),
        S * (1 - A) * g0_scores(theta[g0_pos]),
        S * A * g1_scores(theta[g1_pos]),
        (1 - A) * h_score_weights * h_scores(theta[h_pos]),
        (1 - A) * eta_scores(theta[eta_pos]),
        S - q,
        S / q * (A / e * (Y - g1) + g1 - theta[length(theta) - 4]),
        S / q * ((1 - A) / (1 - e) * (Y - g0) + g0 - theta[length(theta) - 3]),
        S / q * ((1 - A) / (1 - e) * (Y - h) + h - theta[length(theta) - 2]),
        ate_trial - theta[length(theta) - 1],
        ate_r_aware - theta[length(theta)]
      )
    }
  }


estimate_r_aware <-
  function(data,
           outcome_formula,
           treatment_formula,
           participant_formula,
           propensity_score,
           use_no_weights_h = FALSE,
           se_bias_correction = TRUE) {
    
    ############################################################################
    # Do m-estimation 
    ############################################################################
    
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
    
    
    models_r_aware <-
      list(
        e = e_model,
        g0 = g0_model,
        g1 = g1_model,
        h = h_model,
        eta = eta_model
      )
    nparms_r_aware <-
      sum(unlist(lapply(models_r_aware, function(x)
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
      estFUN = ee.r_aware,
      data = data,
      root_control = setup_root_control(start = start_val),
      outer_args = list(
        models = models_r_aware,
        propensity_score = propensity_score,
        use_no_weights_h = use_no_weights_h
      ),
      corrections = list(bias_correction_.fay = correction(fay_bias_correction))
    )
    
    ############################################################################
    
    ############
    # Get estimates for ATE
    est_trial <- roots(out)[nparms_r_aware - 1]
    est_r_aware <- roots(out)[nparms_r_aware]
    
    # Compute standard errors and combined estimator
    se_trial <- sqrt(vcov(out)[nparms_r_aware - 1 , nparms_r_aware - 1])
    se_r_aware <- sqrt(vcov(out)[nparms_r_aware, nparms_r_aware])
    cov_r_aware_trial <- vcov(out)[nparms_r_aware, nparms_r_aware - 1]
    lambda <- (se_trial^2 - cov_r_aware_trial) / (se_trial^2 + se_r_aware^2 - 2 * cov_r_aware_trial)
    est_combined <- (1 - lambda) * est_trial + lambda * est_r_aware
    se_combined <- sqrt((se_trial^2 * se_r_aware^2 - cov_r_aware_trial ^ 2) / (se_trial^2 + se_r_aware^2 - 2 * cov_r_aware_trial))
    
    # Compute bias-corrected standard errors
    corr_se_trial <- sqrt(get_corrections(out)$bias_correction_.fay[nparms_r_aware - 1, nparms_r_aware - 1])
    corr_se_r_aware <- sqrt(get_corrections(out)$bias_correction_.fay[nparms_r_aware, nparms_r_aware])
    corr_cov_r_aware_trial <- get_corrections(out)$bias_correction_.fay[nparms_r_aware, nparms_r_aware - 1]
    corr_se_combined <- sqrt((corr_se_trial^2 * corr_se_r_aware^2 - corr_cov_r_aware_trial ^ 2) / (corr_se_trial^2 + corr_se_r_aware^2 - 2 * corr_cov_r_aware_trial))
    
    ############# 
    #Get estimates for control mean outcome
    est_trial.A0 <- roots(out)[nparms_r_aware - 3]
    est_r_aware.A0 <- roots(out)[nparms_r_aware - 2]
    
    # Compute standard errors and combined estimator
    se_trial.A0 <- sqrt(vcov(out)[nparms_r_aware - 3 , nparms_r_aware - 3])
    se_r_aware.A0 <- sqrt(vcov(out)[nparms_r_aware - 2, nparms_r_aware - 2])
    cov_r_aware_trial.A0 <- vcov(out)[nparms_r_aware - 3, nparms_r_aware - 2]
    lambda.A0 <- (se_trial.A0^2 - cov_r_aware_trial.A0) / (se_trial.A0^2 + se_r_aware.A0^2 - 2 * cov_r_aware_trial.A0)
    est_combined.A0 <- (1 - lambda.A0) * est_trial.A0 + lambda.A0 * est_r_aware.A0
    se_combined.A0 <- sqrt((se_trial.A0^2 * se_r_aware.A0^2 - cov_r_aware_trial.A0 ^ 2) / (se_trial.A0^2 + se_r_aware.A0^2 - 2 * cov_r_aware_trial.A0))
    
    # Compute bias-corrected standard errors
    corr_se_trial.A0 <- sqrt(get_corrections(out)$bias_correction_.fay[nparms_r_aware - 3, nparms_r_aware - 3])
    corr_se_r_aware.A0 <- sqrt(get_corrections(out)$bias_correction_.fay[nparms_r_aware - 2, nparms_r_aware - 2])
    corr_cov_r_aware_trial.A0 <- get_corrections(out)$bias_correction_.fay[nparms_r_aware - 3, nparms_r_aware - 2]
    corr_se_combined.A0 <- sqrt((corr_se_trial.A0^2 * corr_se_r_aware.A0^2 - corr_cov_r_aware_trial.A0 ^ 2) / (corr_se_trial.A0^2 + corr_se_r_aware.A0^2 - 2 * corr_cov_r_aware_trial.A0))
    
    return(
      list(
        est_trial = est_trial,
        est_r_aware = est_r_aware,
        est_combined = est_combined,
        lambda = lambda,
        se_trial = se_trial,
        se_r_aware = se_r_aware,
        se_combined = se_combined,
        ci.ll_trial = est_trial - 1.96 * se_trial,
        ci.ul_trial = est_trial + 1.96 * se_trial,
        ci.ll_r_aware = est_r_aware - 1.96 * se_r_aware,
        ci.ul_r_aware = est_r_aware + 1.96 * se_r_aware,
        ci.ll_combined = est_combined - 1.96 * se_combined,
        ci.ul_combined = est_combined + 1.96 * se_combined,
        corr_se_trial = corr_se_trial,
        corr_se_r_aware = corr_se_r_aware,
        corr_se_combined = corr_se_combined,
        corr_ci.ll_trial = est_trial - 1.96 * corr_se_trial,
        corr_ci.ul_trial = est_trial + 1.96 * corr_se_trial,
        corr_ci.ll_r_aware = est_r_aware - 1.96 * corr_se_r_aware,
        corr_ci.ul_r_aware = est_r_aware + 1.96 * corr_se_r_aware,
        corr_ci.ll_combined = est_combined - 1.96 * corr_se_combined,
        corr_ci.ul_combined = est_combined + 1.96 * corr_se_combined,
        est_trial.A0 = est_trial.A0,
        est_r_aware.A0 = est_r_aware.A0,
        est_combined.A0 = est_combined.A0,
        lambda.A0 = lambda.A0,
        se_trial.A0 = se_trial.A0,
        se_r_aware.A0 = se_r_aware.A0,
        se_combined.A0 = se_combined.A0,
        ci.ll_trial.A0 = est_trial.A0 - 1.96 * se_trial.A0,
        ci.ul_trial.A0 = est_trial.A0 + 1.96 * se_trial.A0,
        ci.ll_r_aware.A0 = est_r_aware.A0 - 1.96 * se_r_aware.A0,
        ci.ul_r_aware.A0 = est_r_aware.A0 + 1.96 * se_r_aware.A0,
        ci.ll_combined.A0 = est_combined.A0 - 1.96 * se_combined.A0,
        ci.ul_combined.A0 = est_combined.A0 + 1.96 * se_combined.A0,
        corr_se_trial.A0 = corr_se_trial.A0,
        corr_se_r_aware.A0 = corr_se_r_aware.A0,
        corr_se_combined.A0 = corr_se_combined.A0,
        corr_ci.ll_trial.A0 = est_trial.A0 - 1.96 * corr_se_trial.A0,
        corr_ci.ul_trial.A0 = est_trial.A0 + 1.96 * corr_se_trial.A0,
        corr_ci.ll_r_aware.A0 = est_r_aware.A0 - 1.96 * corr_se_r_aware.A0,
        corr_ci.ul_r_aware.A0 = est_r_aware.A0 + 1.96 * corr_se_r_aware.A0,
        corr_ci.ll_combined.A0 = est_combined.A0 - 1.96 * corr_se_combined.A0,
        corr_ci.ul_combined.A0 = est_combined.A0 + 1.96 * corr_se_combined.A0
      )
    )
    
  }
