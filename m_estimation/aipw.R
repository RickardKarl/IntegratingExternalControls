ee.trial <- function(data, models, propensity_score = NULL) {
  # Get design matrix
  Xe  <-
    grab_design_matrix(data, rhs_formula = grab_fixed_formula(models$e))
  Xg0 <-
    grab_design_matrix(data, rhs_formula = grab_fixed_formula(models$g0))
  Xg1 <-
    grab_design_matrix(data, rhs_formula = grab_fixed_formula(models$g1))
  
  # Get indices for parameters
  e_pos  <- 1:ncol(Xe)
  g0_pos <- (max(e_pos)  + 1):(max(e_pos)  + ncol(Xg0))
  g1_pos <- (max(g0_pos) + 1):(max(g0_pos) + ncol(Xg1))
  
  # Grab estimating function from each model object
  e_scores  <- grab_psiFUN(models$e, data)
  g0_scores <- grab_psiFUN(models$g0, data)
  g1_scores <- grab_psiFUN(models$g1, data)
  
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
    q  <- theta[length(theta) - 3]
    ate_trial <-
      (theta[length(theta) - 2] - theta[length(theta) - 1])
    c(
      S * e_scores(theta[e_pos]),
      S * (1 - A) * g0_scores(theta[g0_pos]),
      S * A * g1_scores(theta[g1_pos]),
      S - q,
      S / q * (A / e * (Y - g1) + g1 - theta[length(theta) - 2]),
      S / q * ((1 - A) / (1 - e) * (Y - g0) + g0 - theta[length(theta) - 1]),
      S * (ate_trial - theta[length(theta)])
    )
  }
}

estimate_aipw <-
  function(data,
           outcome_formul,
           treatment_formula,
           propensity_score) {
    
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
    
    g0_model <-
      glm(outcome_formula,
          subset = (A == 0 & S == 1),
          data = data)
    
    g1_model <-
      glm(outcome_formula,
          subset = (A == 1 & S == 1),
          data = data)
    
    models <-
      list(e = e_model,
           g0 = g0_model,
           g1 = g1_model)
    
    
    start_val <-
      c(
        coef(e_model),
        coef(g0_model),
        coef(g1_model),
        mean(data$S),
        mean(data$Y[(data$A == 1) & (data$S == 1)]),
        mean(data$Y[(data$A == 0) & (data$S == 1)]),
        mean(data$Y[(data$A == 1) &
                      (data$S == 1)]) -  mean(data$Y[(data$A == 0) &
                                                       (data$S == 1)])
      )
    nparms <- length(start_val)
    
    out <- m_estimate(
      estFUN = ee.trial,
      data = data,
      root_control = setup_root_control(start = start_val),
      outer_args = list(models = models, propensity_score = propensity_score),
      corrections = list(bias_correction_.fay = correction(fay_bias_correction))
    )
    
    ############################################################################
  
    # get estimates of ATE
    estimate = roots(out)[nparms]
    se = sqrt(vcov(out)[nparms, nparms])
    corr_se = sqrt(get_corrections(out)$bias_correction_.fay[nparms, nparms])
    
    # get estimates of control mean outcome
    estimate.A0 = roots(out)[nparms - 1]
    se.A0 = sqrt(vcov(out)[nparms - 1, nparms - 1])
    corr_se.A0 =  sqrt(get_corrections(out)$bias_correction_.fay[nparms - 1, nparms - 1])
    
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