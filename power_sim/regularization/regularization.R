
library(glmnet)


truncated_lines <- function(x, kappa){
  (x - kappa)*(x > kappa)
}


S_lasso_ridge <- function(data, significance_level = 0.05) {
  
  States = unique(data$State)
  data[["id"]] <- paste0(data$Rep, data$Charge, data$Experimental_state)
  Test = aic = loglik = Test_statistic = p_value = rep(NA, 2)
  
  knots <- unique(setdiff(data$Exposure, c(max(data$Exposure), min(data$Exposure))))
  X <- sapply(knots, function(k) {
    truncated_lines(data$Exposure, k)
  })
  colnames(X) <- c(paste0("knot_", as.character(knots)))
  
  
  #LASSO
  
  cv_fit <- glmnet(X, data[["Mass"]], alpha = 1, lambda = 0.001)
  coefs <- coefficients(cv_fit)
  X_reduced <- cbind(intercept = 1, X)[, which(coefs != 0)]
  
  # continuous, identity
  model = lmerTest::lmer(Mass ~ Exposure*Experimental_state + (1|id) + X_reduced,
                         data = data,
                         REML = FALSE)
  model_reduced = lmerTest::lmer(Mass ~ Exposure + (1|id) + X_reduced,
                                 data = data,
                                 REML = FALSE)
  result = anova(model, model_reduced)
  aic[1] = AIC(model)
  loglik[1] = as.numeric(logLik(model))
  Test_statistic[1] = result$Chisq[2]
  p_value[1] = result$`Pr(>Chisq)`[2]
  
  Test[1] = "LASSO_knots_random_intercept_id"
  
  
  #RIDGE
  
  cv_fit <- glmnet(X, data[["Mass"]], alpha = 0, lambda = 0.001)
  coefs <- coefficients(cv_fit)
  X_reduced <- cbind(intercept = 1, X)[, which(abs(coefs) >= 2*10^(-5))]
  
  # continuous, identity
  model = lmerTest::lmer(Mass ~ Exposure*Experimental_state + (1|id) + X_reduced,
                         data = data,
                         REML = FALSE)
  model_reduced = lmerTest::lmer(Mass ~ Exposure + (1|id) + X_reduced,
                                 data = data,
                                 REML = FALSE)
  result = anova(model, model_reduced)
  aic[2] = AIC(model)
  loglik[2] = as.numeric(logLik(model))
  Test_statistic[2] = result$Chisq[2]
  p_value[2] = result$`Pr(>Chisq)`[2]
  
  Test[2] = "RIDGE_knots_random_intercept_id"
  
  data.table::data.table(Test = Test,
                         State_1 = States[1],
                         State_2 = States[2],
                         Test_statistic = Test_statistic,
                         P_value = p_value,
                         Significant_difference = (p_value <= significance_level),
                         Time = NA,
                         Transformation = "identity",
                         AIC = aic,
                         logLik = loglik)
}

