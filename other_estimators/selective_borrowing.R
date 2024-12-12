# need to install package with "devtools::install_github("Gaochenyin/SelectiveIntegrative")"

# this is a wrapper for their method from
# Integrating Randomized Placebo-Controlled Trial Data with External Controls:
# A Semiparametric Approach with Selective Borrowing


library(SelectiveIntegrative)


selective_borrowing.ate <-
  function(data,
           outcome_formula) {
    # Get design matrix of outcome_formula without intercept
    
    
    
    Y <- data$Y
    X <- data.frame(model.matrix(outcome_formula, data)[, -1])
    A <- data$A
    S <- data$S
    
    
    
    data_rt <-
      list(Y = Y[S == 1], A = A[S == 1], X = X[S ==
                                                 1,])
    data_ec <-
      list(Y = Y[S == 0], A = A[S == 0], X = X[S ==
                                                 0,])
    {
      sink("/dev/null") # to silence print statement inside srEC
      
      res <-      srEC(
        data_rt = data_rt,
        data_ec = list(data_ec),
        rt.ctrl = caret::trainControl(method = 'cv', number = 3),
        hc.ctrl = caret::trainControl(method = 'cv', number = 3),
        method = 'glm'
      )
      sink()
      }
    
    # get the selective integrative estimator components
    est <- res$est$ACW.final
    se <- res$sd$ACW.final / sqrt(res$n_c)
    return(data.frame(
      est = est,
      se = se,
      ci.ll = est - 1.96 * se,
      ci.ul = est + 1.96 * se,
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
    ))
  }
