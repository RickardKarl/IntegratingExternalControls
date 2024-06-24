library(progress) # progress bar
library(geex)

source('m_estimation/ipw.R')
source('m_estimation/randomization-aware.R')
source('m_estimation/pooling.R')

source('m_estimation/test-then-pool.R')
source('dynamic_borrowing/mixture_prior.R')
source('dynamic_borrowing/selective_borrowing.R')


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
  # If we do not provide arguments through command-line
  in_n1 <- 100
  in_n0 <- 100
  in_delta <- 0
  in_corr_spec <- 1
  seed = 20240624
}

set.seed(seed)
iterations <- 250

dataset = 'data/polynomial_dgp.R'
source(dataset)
tparam <- true_ate() # true ATE

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
    'robin' = function(data)
      estimate_robin(
        data,
        outcome_formula = outcome_formula,
        treatment_formula = treatment_formula,
        participant_formula = participant_formula
      ),
    'ipw' = function(data)
      estimate_ipw(data),
    'pooling' = function(data)
      estimate_pooling(
        data,
        outcome_formula = outcome_formula,
        treatment_formula = treatment_formula,
        participant_formula = participant_formula
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
        participant_formula = participant_formula
      )
  )

method_list <-
  available_estimators[c(
    'robin',
    'ipw',
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
    'coverage'
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
  
  for (j in 1:length(method_list)) {
    out <- method_list[[j]](data)
    
    
    if (names(method_list)[j] == 'robin') {
      estimates = matrix(nrow = 3)
      estimates[1] <- out$est_trial
      estimates[2] <- out$est_robin
      estimates[3] <- out$est_combined
      
      se = matrix(nrow = 3)
      se[1] <- out$se_trial
      se[2] <- out$se_robin
      se[3] <- out$se_combined
      
      ci_ll = matrix(nrow = 3)
      ci_ul = matrix(nrow = 3)
      ci_ll[1] <- out$ci.ll_trial
      ci_ul[1] <- out$ci.ul_trial
      ci_ll[2] <- out$ci.ll_robin
      ci_ul[2] <- out$ci.ul_robin
      ci_ll[3] <- out$ci.ll_combined
      ci_ul[3] <- out$ci.ul_combined
      
      coverage = matrix(nrow = 3)
      coverage[1] <- (ci_ll[1]  < tparam) & (tparam <  ci_ul[1])
      coverage[2] <- (ci_ll[2]  < tparam) & (tparam <  ci_ul[2])
      coverage[3] <- (ci_ll[3]  < tparam) & (tparam <  ci_ul[3])
      
      
      
      est_names <- c('trial', 'robin', 'combined')
      for (j in 1:3) {
        res_df[nrow(res_df) + 1,] <-
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
            coverage[j]
          )
      }
    } else {
      res_df[nrow(res_df) + 1,] <-
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
          (out$ci.ll  < tparam) & (tparam <  out$ci.ul)
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
