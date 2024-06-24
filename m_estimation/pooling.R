#source('m_estimation/bias_correction.R')

ee.pooling <- function(data, models) {
  S <- data$S
  A <- data$A
  Y <- data$Y
  
  # Get design matrix
  Xe  <-
    grab_design_matrix(data, rhs_formula = grab_fixed_formula(models$e))
  Xg0 <-
    grab_design_matrix(data, rhs_formula = grab_fixed_formula(models$g0))
  Xg1 <-
    grab_design_matrix(data, rhs_formula = grab_fixed_formula(models$g1))
  Xeta <-
    grab_design_matrix(data, rhs_formula = grab_fixed_formula(models$eta))
  
  # Get indices for parameters
  e_pos  <- 1:ncol(Xe)
  g0_pos <-  (max(e_pos)  + 1):(max(e_pos)  + ncol(Xg0))
  g1_pos <-  (max(g0_pos) + 1):(max(g0_pos) + ncol(Xg1))
  eta_pos <- (max(g1_pos)  + 1):(max(g1_pos)  + ncol(Xeta))
  
  # Grab estimating function from each model object
  e_scores  <- grab_psiFUN(models$e, data)
  g0_scores <- grab_psiFUN(models$g0, data)
  g1_scores <- grab_psiFUN(models$g1, data)
  eta_scores <- grab_psiFUN(models$eta, data)
  
  function(theta) {
    e  <- plogis(Xe %*% theta[e_pos])[, 1]
    g0 <- Xg0 %*% theta[g0_pos]
    g1 <- Xg1 %*% theta[g1_pos]
    eta <- plogis(Xeta %*% theta[eta_pos])[, 1]
    q  <- theta[length(theta) - 3]
    weights.A0 <-
      (S * (1 - A) + (1 - S)) / (eta * (1 - e) + 1 - eta) * eta
    ate_pooling <-
      theta[length(theta) - 2] - theta[length(theta) - 1]
    c(
      S * e_scores(theta[e_pos]),
      (1 - A) * g0_scores(theta[g0_pos]),
      S * A * g1_scores(theta[g1_pos]),
      eta_scores(theta[eta_pos]),
      S - q,
      S / q * ((A / e * (Y - g1) + g1) - theta[length(theta) - 2]),
      1 / q * (weights.A0 * (Y - g0) + S * (g0 - theta[length(theta) - 1])),
      ate_pooling - theta[length(theta)]
    )
    
  }
}

estimate_pooling <-
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
      glm(participant_formula,
          data = data,
          family = binomial)
    
    g0_model <-
      glm(outcome_formula,
          subset = (A == 0),
          data = data)
    
    g1_model <-
      glm(outcome_formula,
          subset = (A == 1 & S == 1),
          data = data)
    
    models <-
      list(e = e_model,
           g0 = g0_model,
           g1 = g1_model,
           eta = eta_model)
    nparms <-
      sum(unlist(lapply(models, function(x)
        length(coef(
          x
        ))))) + 4
    
    start_val <-
      c(
        coef(e_model),
        coef(g0_model),
        coef(g1_model),
        coef(eta_model),
        mean(data$S),
        mean(data$Y[(data$A == 1) & (data$S == 1)]),
        mean(data$Y[(data$A == 0) & (data$S == 1)]),
        mean(data$Y[(data$A == 1) &
                      (data$S == 1)]) -  mean(data$Y[(data$A == 0) &
                                                       (data$S == 1)])
      )
    
    out <- m_estimate(
      estFUN = ee.pooling,
      data = data,
      root_control = setup_root_control(start = start_val),
      outer_args = list(models = models),
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
