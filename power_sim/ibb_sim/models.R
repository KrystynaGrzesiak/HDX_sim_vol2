library(lmerTest)
library(glmnet)



truncated_lines <- function(x, knots){
  sapply(knots, function(kappa) {
    (x - kappa)*(x > kappa)
  })
}


test_Eq1 <- function(data, significance_level = 0.05) {
  
  n_knots <- 1:5
  
  #parameters
  States <- unique(data[["State"]])
  aic <- loglik <- Test_statistic <- p_value <- rep(NA, length(n_knots))
  Times <- unique(data[["Exposure"]])
  
  #semiparametric
  t <- (Times[1:(length(Times)-1)] + Times[2:length(Times)])/2
  knots <- unique(setdiff(t, c(max(t), min(t))))
  X <- truncated_lines(data[["Exposure"]], knots)
  colnames(X) <- c(paste0("knot_", as.character(knots)))
  cv_fit <- glmnet(X, data[["Mass"]], alpha = 0, lambda = 0.001)
  coefs <- as.numeric(coefficients(cv_fit))
  
  #model
  do.call(rbind, lapply(n_knots, function(i) {
    
    #choice of knots
    X_reduced <- cbind(intercept = 1, X)[, abs(coefs) %in% tail(sort(abs(coefs)), i)]
    
    suppressMessages({
      model <- lmer(Mass ~ Exposure*State + (1|Exposure) + (1|State) + X_reduced,
                    data = data,
                    REML = FALSE)
      model_reduced <- lmer(Mass ~ Exposure + (1|Exposure) + X_reduced,
                            data = data,
                            REML = FALSE)
    })
    
    result <- anova(model, model_reduced)
    aic <- AIC(model)
    loglik <- as.numeric(logLik(model))
    Test_statistic <- result$Chisq[2]
    p_value <- result$`Pr(>Chisq)`[2]
    
    data.table(Test = paste0("Eq1_knots_", i),
               State_1 = States[1],
               State_2 = States[2],
               Test_statistic = Test_statistic,
               P_value = p_value,
               Significant_difference = (p_value <= significance_level),
               Time = NA,
               Transformation = "identity",
               AIC = aic,
               logLik = loglik)
  }))
}


test_Eq2 <- function(data, significance_level = 0.05) {
  
  n_knots <- 1:5
  
  #parameters
  States <- unique(data[["State"]])
  aic <- loglik <- Test_statistic <- p_value <- rep(NA, length(n_knots))
  Times <- unique(data[["Exposure"]])
  
  #semiparametric
  t <- (Times[1:(length(Times)-1)] + Times[2:length(Times)])/2
  knots <- unique(setdiff(t, c(max(t), min(t))))
  X <- truncated_lines(data[["Exposure"]], knots)
  colnames(X) <- c(paste0("knot_", as.character(knots)))
  cv_fit <- glmnet(X, data[["Mass"]], alpha = 0, lambda = 0.001)
  coefs <- as.numeric(coefficients(cv_fit))
  
  #model
  do.call(rbind, lapply(n_knots, function(i) {
    
    #choice of knots
    X_reduced <- cbind(intercept = 1, X)[, abs(coefs) %in% tail(sort(abs(coefs)), i)]
    
    suppressMessages({
      model <- lmer(Mass ~ Exposure*State + (1|Exposure) + X_reduced,
                    data = data,
                    REML = FALSE)
      model_reduced <- lmer(Mass ~ Exposure + (1|Exposure) + X_reduced,
                            data = data,
                            REML = FALSE)
    })
    
    result <- anova(model, model_reduced)
    aic <- AIC(model)
    loglik <- as.numeric(logLik(model))
    Test_statistic <- result$Chisq[2]
    p_value <- result$`Pr(>Chisq)`[2]
    
    data.table(Test = paste0("Eq2_knots_", i),
               State_1 = States[1],
               State_2 = States[2],
               Test_statistic = Test_statistic,
               P_value = p_value,
               Significant_difference = (p_value <= significance_level),
               Time = NA,
               Transformation = "identity",
               AIC = aic,
               logLik = loglik)
  }))
}



test_Eq3 <- function(data, significance_level = 0.05) {
  
  n_knots <- 1:5
  
  #parameters
  States <- unique(data[["State"]])
  aic <- loglik <- Test_statistic <- p_value <- rep(NA, length(n_knots))
  Times <- unique(data[["Exposure"]])
  
  #semiparametric
  t <- (Times[1:(length(Times)-1)] + Times[2:length(Times)])/2
  knots <- unique(setdiff(t, c(max(t), min(t))))
  X <- truncated_lines(data[["Exposure"]], knots)
  colnames(X) <- c(paste0("knot_", as.character(knots)))
  cv_fit <- glmnet(X, data[["Mass"]], alpha = 0, lambda = 0.001)
  coefs <- as.numeric(coefficients(cv_fit))
  
  #model
  do.call(rbind, lapply(n_knots, function(i) {
    
    #choice of knots
    X_reduced <- cbind(intercept = 1, X)[, abs(coefs) %in% tail(sort(abs(coefs)), i)]
    
    suppressMessages({
      model <- lmer(Mass ~ Exposure*State + (1|Exposure) + X_reduced + (Exposure || State),
                    data = data,
                    REML = FALSE)
      model_reduced <- lmer(Mass ~ Exposure + (1|Exposure) + X_reduced,
                            data = data,
                            REML = FALSE)
    })
    
    result <- anova(model, model_reduced)
    aic <- AIC(model)
    loglik <- as.numeric(logLik(model))
    Test_statistic <- result$Chisq[2]
    p_value <- result$`Pr(>Chisq)`[2]
    
    data.table(Test = paste0("Eq2_knots_", i),
               State_1 = States[1],
               State_2 = States[2],
               Test_statistic = Test_statistic,
               P_value = p_value,
               Significant_difference = (p_value <= significance_level),
               Time = NA,
               Transformation = "identity",
               AIC = aic,
               logLik = loglik)
  }))
}



tests_fun <- list(test_Eq1, test_Eq2, test_Eq3)



