
library(glmnet)
library(data.table)

truncated_lines <- function(x, knots){
  sapply(knots, function(kappa) {
    (x - kappa)*(x > kappa)
  })
}

test_semiparametric <- function(data, significance_level = 0.05) {
  
  States <- unique(data[["State"]])
  data[["id"]] <- paste0(data[["Rep"]], data[["Charge"]], data[["Experimental_state"]])
  Test <- aic <- loglik <- Test_statistic <- p_value <- NA
  
  Times <- unique(data[["Exposure"]])
  
  if(length(Times) > 2) {
    
    knots <- unique(setdiff(data[["Exposure"]], c(max(data[["Exposure"]]), min(data[["Exposure"]]))))
    X <- truncated_lines(data[["Exposure"]], knots)

    colnames(X) <- c(paste0("knot_", as.character(knots)))
    
    cv_fit <- glmnet(X, data[["Mass"]], alpha = 0, lambda = 0.001)
    coefs <- coefficients(cv_fit)
    X_reduced <- cbind(intercept = 1, X)[, which(as.logical(abs(coefs) >= 2*10^(-5)))]
    
  } else {
    if(length(Times) == 2) {
      
      X_reduced <- rep(0, nrow(data))
      
    }else{
      
      stop("You must provide more than one time point.")
    }
  }
  
  suppressMessages({
    model <- lmer(Mass ~ Exposure*State  + (1|Exposure) + X_reduced,
                  data = data,
                  REML = FALSE)
    model_reduced <- lmer(Mass ~ Exposure  + (1|Exposure) + X_reduced,
                          data = data,
                          REML = FALSE)
  })
  
  result <- anova(model, model_reduced)
  aic <- AIC(model)
  loglik <- as.numeric(logLik(model))
  Test_statistic <- result$Chisq[2]
  p_value <- result$`Pr(>Chisq)`[2]
  
  Test <- "Semiparametric model"
  
  data.table(Test = Test,
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
