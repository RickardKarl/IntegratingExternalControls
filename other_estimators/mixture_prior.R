check_diagnostics <- function(mcmc_model, n_chains) {
  diagnostics <- hb_convergence(mcmc_model)
  if (diagnostics$max_rhat >= 1.01) {
    print('max_rhat is too large')
    return(FALSE)
  } else if (diagnostics$min_ess_bulk < 100 * n_chains |
             diagnostics$min_ess_tail < 100 * n_chains) {
    print("min_ess_bulk or min_ess_tail is too small")
    return(FALSE)
  } else {
    return(TRUE)
  }
}

mixture_prior.ate <- function(data,
                              covariates,
                              n_chains = 4,
                              n_adapt = 2e3,
                              n_warmup = 2e3,
                              n_iterations = 4e3) {
  require(historicalborrow)
  

  data$id = rownames(data)

  
  # Get hyperparameters
  hyperparameters <- hb_mcmc_mixture_hyperparameters(
    data = data,
    response = 'Y',
    study = 'S',
    study_reference = 1,
    group = 'A',
    group_reference = 0,
    patient = 'id',
  )
  
  # Get posterior with mixture prior
  data_mixture <- dplyr::filter(data, data$S == 1)
  mcmc_mixture <- hb_mcmc_mixture(
    data = data_mixture,
    # only analyze current study
    response = 'Y',
    study = 'S',
    study_reference = 1,
    group = 'A',
    group_reference = 0,
    patient = 'id',
    covariates = covariates,
    m_omega = hyperparameters$m_omega,
    s_omega = hyperparameters$s_omega,
    p_omega = rep(1 / nrow(hyperparameters), nrow(hyperparameters)),
    n_chains = n_chains,
    n_adapt = n_adapt,
    n_warmup = n_warmup,
    n_iterations = n_iterations
  )
  
  
  # Test diagnostics (this should return true)
  if (!check_diagnostics(mcmc_mixture, n_chains)) {
    warning('MCMC went wrong!')
  }
  
  #### Get summaries
  
  summary_mixture <- hb_summary(
    mcmc = mcmc_mixture,
    data = data_mixture,
    response = 'Y',
    study = 'S',
    study_reference = 1,
    group = 'A',
    group_reference = 0,
    patient = 'id',
    covariates = covariates,
  )
  #### Get summary statistics
  return(
    data.frame(
      est = summary_mixture$diff_mean[2],
      se = NaN,
      ci.ll = summary_mixture$diff_lower[2],
      ci.ul = summary_mixture$diff_upper[2],
      corr_se = NaN,
      corr_ci.ll = NaN,
      corr_ci.ul = NaN,
      est.A0 = NaN,
      se.A0 = NaN,
      ci.ll.A0 = NaN,
      ci.ul.A0 = NaN,
      corr_se.A0 = NaN,
      corr_ci.ll.A0 = NaN,
      corr_ci.ul.A0 = NaN
    )
  )
}
