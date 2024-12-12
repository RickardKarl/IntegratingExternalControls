library(progress) # progress bar
library(geex)

source('m_estimation/ipw.R')
source('m_estimation/r_aware.R')
source('m_estimation/pooling.R')

source('m_estimation/test-then-pool.R')
source('other_estimators/mixture_prior.R')
source('other_estimators/selective_borrowing.R')


ts <- format(Sys.time(), "%d%m_%H%M%S")
print(paste('Timestamp for experiment:', ts))

########################
# CONFIG FOR EXPERIMENT
########################

# Get input arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) >= 5) {
  in_n1 <- as.integer(args[1])
  in_n0 <- as.integer(args[2])
  in_delta <- as.numeric(args[3])
  in_corr_spec <- as.integer(args[4])
  seed <- as.integer(args[5])
  
  cat('n1', in_n1, '\n')
  cat('n0', in_n0, '\n')
  cat('delta', in_delta, '\n')
  cat('correct model spec.', in_corr_spec, '\n')
  cat('seed', seed, '\n')
} else {
  in_n1 <- 50
  in_n0 <- 200
  in_delta <- 0
  in_corr_spec <- 1
  seed = 50
}

set.seed(seed)
iterations <- 250

dataset = 'data/polynomial_dgp.R'
source(dataset)
tparam <- true_ate() # true ATE
tparam.A0 <- true_control_mean()
propensity_score = 1 / 2

##################
# SETUP METHODS
##################



if (in_corr_spec == 1) {
  outcome_formula <- formula(
    Y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + I(X1 ^ 2) + I(X2 ^ 2) + I(X3 ^ 2) + I(X4 ^ 2) + I(X5 ^ 2) + I(X6 ^ 2) + I(X7 ^ 2) + I(X8 ^ 2) + I(X9 ^ 2) + I(X10 ^ 2)
  )
  treatment_formula <-
    formula(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10)
  participant_formula <-
    formula(S ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10)
  covariates <-
    c('X1', 'X2', 'X3', 'X4', 'X5', 'X6', 'X7', 'X8', 'X9', 'X10')
} else {
  outcome_formula <-
    formula(Y ~ X1 + X2 + X3 + X4 + X5)
  treatment_formula <- formula(A ~ X1 + X2 + X3 + X4 + X5)
  participant_formula <- formula(S ~ X1 + X2 + X3 + X4 + X5)
  covariates <- c('X1', 'X2', 'X3', 'X4', 'X5')
}


available_estimators <-
  c(
    'standard' = function(data)
      estimate_r_aware(
        data,
        outcome_formula = outcome_formula,
        treatment_formula = treatment_formula,
        participant_formula = participant_formula,
        propensity_score = propensity_score
      ),
    'unweighted' = function(data)
      estimate_r_aware(
        data,
        outcome_formula = outcome_formula,
        treatment_formula = treatment_formula,
        participant_formula = participant_formula,
        propensity_score = propensity_score,
        use_no_weights_h = TRUE
      ),
    'estimated_propensity' = function(data)
      estimate_r_aware(
        data,
        outcome_formula = outcome_formula,
        treatment_formula = treatment_formula,
        participant_formula = participant_formula,
        propensity_score = NULL
      ),
    'ipw' = function(data)
      estimate_ipw(data, propensity_score = propensity_score),
    'pooling' = function(data)
      estimate_pooling(
        data,
        outcome_formula = outcome_formula,
        treatment_formula = treatment_formula,
        participant_formula = participant_formula,
        propensity_score = propensity_score
      ),
    'dynamic_borrowing' = function(data)
      mixture_prior.ate(data, covariates),
    'selective_borrowing' = function(data)
      selective_borrowing.ate(data, outcome_formula = outcome_formula),
    'test-then-pool' = function(data)
      test_then_pool(
        data,
        outcome_formula = outcome_formula,
        treatment_formula = treatment_formula,
        participant_formula = participant_formula,
        propensity_score = propensity_score
      )
  )

r_aware_method <- c('standard', 'unweighted', 'estimated_propensity')

method_list <-
  available_estimators[c(
    'ipw',
    'standard',
    'unweighted',
    'estimated_propensity',
    'pooling',
    'test-then-pool',
    'selective_borrowing',
    'dynamic_borrowing'
  )]

##################
# START EXPERIMENT
##################

columns <-
  c(
    'iterations',
    'seed',
    'delta',
    'n1',
    'n0',
    'modelspec',
    'dataset',
    'estimator',
    'estimate',
    'se',
    'ci_ll',
    'ci_ul',
    'coverage',
    'corr_se',
    'corr_ci_ll',
    'corr_ci_ul',
    'corr_coverage',
    'lambda',
    'estimate.A0',
    'se.A0',
    'ci_ll.A0',
    'ci_ul.A0',
    'coverage.A0',
    'corr_se.A0',
    'corr_ci_ll.A0',
    'corr_ci_ul.A0',
    'corr_coverage.A0'
  )
res_df <- data.frame(matrix(ncol = length(columns), nrow = 0))
colnames(res_df) <- columns

# Initializes the progress bar and time tracking
pb <-
  progress_bar$new(
    format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
    total = iterations,
    complete = "=",
    # Completion bar character
    incomplete = "-",
    # Incomplete bar character
    current = ">",
    # Current bar character
    clear = FALSE,
    # If TRUE, clears the bar when finish
    width = 100
  )      # Width of the progress bar


for (i in 1:iterations)
{
  # Updates the current state
  pb$tick()
  
  data <- dgp(in_n1, in_n0, in_delta)
  
  for (j in 1:length(method_list))
  {
    out <- method_list[[j]](data)
    
    
    if (names(method_list)[j] %in% r_aware_method)
    {
      estimates = matrix(nrow = 3)
      estimates[1] <- out$est_trial
      estimates[2] <- out$est_r_aware
      estimates[3] <- out$est_combined
      
      se = matrix(nrow = 3)
      se[1] <- out$se_trial
      se[2] <- out$se_r_aware
      se[3] <- out$se_combined
      
      corr_se = matrix(nrow = 3)
      corr_se[1] <- out$corr_se_trial
      corr_se[2] <- out$corr_se_r_aware
      corr_se[3] <- out$corr_se_combined
      
      ci_ll = matrix(nrow = 3)
      ci_ul = matrix(nrow = 3)
      ci_ll[1] <- out$ci.ll_trial
      ci_ul[1] <- out$ci.ul_trial
      ci_ll[2] <- out$ci.ll_r_aware
      ci_ul[2] <- out$ci.ul_r_aware
      ci_ll[3] <- out$ci.ll_combined
      ci_ul[3] <- out$ci.ul_combined
      
      corr_ci_ll = matrix(nrow = 3)
      corr_ci_ul = matrix(nrow = 3)
      corr_ci_ll[1] <- out$corr_ci.ll_trial
      corr_ci_ul[1] <- out$corr_ci.ul_trial
      corr_ci_ll[2] <- out$corr_ci.ll_r_aware
      corr_ci_ul[2] <- out$corr_ci.ul_r_aware
      corr_ci_ll[3] <- out$corr_ci.ll_combined
      corr_ci_ul[3] <- out$corr_ci.ul_combined
      
      
      coverage = matrix(nrow = 3)
      coverage[1] <- (ci_ll[1]  < tparam) & (tparam <  ci_ul[1])
      coverage[2] <- (ci_ll[2]  < tparam) & (tparam <  ci_ul[2])
      coverage[3] <- (ci_ll[3]  < tparam) & (tparam <  ci_ul[3])
      
      corr_coverage = matrix(nrow = 3)
      corr_coverage[1] <- (corr_ci_ll[1]  < tparam) & (tparam <  corr_ci_ul[1])
      corr_coverage[2] <- (corr_ci_ll[2]  < tparam) & (tparam <  corr_ci_ul[2])
      corr_coverage[3] <- (corr_ci_ll[3]  < tparam) & (tparam <  corr_ci_ul[3])
      
      
      lambda = matrix(nrow = 3)
      lambda[1] = NaN
      lambda[2] = NaN
      lambda[3] = out$lambda
      
      
      estimates.A0 = matrix(nrow = 3)
      estimates.A0[1] <- out$est_trial.A0
      estimates.A0[2] <- out$est_r_aware.A0
      estimates.A0[3] <- out$est_combined.A0
      
      se.A0 = matrix(nrow = 3)
      se.A0[1] <- out$se_trial.A0
      se.A0[2] <- out$se_r_aware.A0
      se.A0[3] <- out$se_combined.A0
      
      corr_se.A0 = matrix(nrow = 3)
      corr_se.A0[1] <- out$corr_se_trial.A0
      corr_se.A0[2] <- out$corr_se_r_aware.A0
      corr_se.A0[3] <- out$corr_se_combined.A0
      
      
      ci_ll.A0 = matrix(nrow = 3)
      ci_ul.A0 = matrix(nrow = 3)
      ci_ll.A0[1] <- out$ci.ll_trial.A0
      ci_ul.A0[1] <- out$ci.ul_trial.A0
      ci_ll.A0[2] <- out$ci.ll_r_aware.A0
      ci_ul.A0[2] <- out$ci.ul_r_aware.A0
      ci_ll.A0[3] <- out$ci.ll_combined.A0
      ci_ul.A0[3] <- out$ci.ul_combined.A0
      
      corr_ci_ll.A0 = matrix(nrow = 3)
      corr_ci_ul.A0 = matrix(nrow = 3)
      corr_ci_ll.A0[1] <- out$corr_ci.ll_trial.A0
      corr_ci_ul.A0[1] <- out$corr_ci.ul_trial.A0
      corr_ci_ll.A0[2] <- out$corr_ci.ll_r_aware.A0
      corr_ci_ul.A0[2] <- out$corr_ci.ul_r_aware.A0
      corr_ci_ll.A0[3] <- out$corr_ci.ll_combined.A0
      corr_ci_ul.A0[3] <- out$corr_ci.ul_combined.A0
      
      coverage.A0 = matrix(nrow = 3)
      coverage.A0[1] <- (ci_ll.A0[1]  < tparam.A0) & (tparam.A0 <  ci_ul.A0[1])
      coverage.A0[2] <- (ci_ll.A0[2]  < tparam.A0) & (tparam.A0 <  ci_ul.A0[2])
      coverage.A0[3] <- (ci_ll.A0[3]  < tparam.A0) & (tparam.A0 <  ci_ul.A0[3])
      
      corr_coverage.A0 = matrix(nrow = 3)
      corr_coverage.A0[1] <- (corr_ci_ll.A0[1]  < tparam.A0) & (tparam.A0 <  corr_ci_ul.A0[1])
      corr_coverage.A0[2] <- (corr_ci_ll.A0[2]  < tparam.A0) & (tparam.A0 <  corr_ci_ul.A0[2])
      corr_coverage.A0[3] <- (corr_ci_ll.A0[3]  < tparam.A0) & (tparam.A0 <  corr_ci_ul.A0[3])
      
      
      tag <- names(method_list)[j]
      est_names <- c(
        paste('trial', tag, sep  =  '.')
        ,
        paste('r_aware', tag, sep  =  '.')
        ,
        paste('combined', tag, sep  =  '.')
      )
      for (j in 1:3)
      {
        res_df[nrow(res_df) + 1, ] <-
          c(
            iterations,
            seed,
            in_delta,
            in_n1,
            in_n0,
            in_corr_spec,
            dataset,
            est_names[j],
            estimates[j],
            se[j],
            ci_ll[j],
            ci_ul[j],
            coverage[j],
            corr_se[j],
            corr_ci_ll[j],
            corr_ci_ul[j],
            corr_coverage[j],
            lambda[j],
            estimates.A0[j],
            se.A0[j],
            ci_ll.A0[j],
            ci_ul.A0[j],
            coverage.A0[j],
            corr_se.A0[j],
            corr_ci_ll.A0[j],
            corr_ci_ul.A0[j],
            corr_coverage.A0[j]
          )
      }
    } else
    {
      res_df[nrow(res_df) + 1, ] <-
        c(
          iterations,
          seed,
          in_delta,
          in_n1,
          in_n0,
          in_corr_spec,
          dataset,
          names(method_list)[j],
          out$est,
          out$se,
          out$ci.ll,
          out$ci.ul,
          (out$ci.ll  < tparam) & (tparam <  out$ci.ul),
          out$corr_se,
          out$corr_ci.ll,
          out$corr_ci.ul,
          (out$corr_ci.ll  < tparam) & (tparam <  out$corr_ci.ul),
          NaN,
          out$est.A0,
          out$se.A0,
          out$ci.ll.A0,
          out$ci.ul.A0,
          (out$ci.ll.A0  < tparam.A0) &
            (tparam.A0 <  out$ci.ul.A0),
          out$corr_se.A0,
          out$corr_ci.ll.A0,
          out$corr_ci.ul.A0,
          (out$corr_ci.ll.A0  < tparam.A0) & (tparam.A0 <  out$corr_ci.ul.A0)
        )
    }
  }
  
  # Save current iterations to file
  write.csv(
    res_df,
    paste('output_folder/',
          ts,
          '-',
          seed,
          '.csv',
          sep = ''),
    row.names = FALSE,
    quote = FALSE
  )
}
