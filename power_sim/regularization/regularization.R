
library(glmnet)


truncated_lines <- function(x, kappa){
  (x - kappa)*(x > kappa)
}


### ridge

S_reg <- function(data, significance_level = 0.05) {
  
  States = unique(data$State)
  
  aic = loglik = Test_statistic = p_value = rep(NA, 2)
  data[["Exposure"]] = log(data[["Exposure"]] + 1)
  
  knots <- unique(setdiff(data$Exposure, c(max(data$Exposure), min(data$Exposure))))
  
  X <- sapply(1:length(knots), function(knot) {
    truncated_lines(data$Exposure, knots[knot])
  })
  
  design_matrix <- cbind(data[["Experimental_state"]] == "A", 
                         X, 
                         data[["Exposure"]]*(data[["Experimental_state"]] == "A"))
  
  colnames(design_matrix) <- c("State", as.character(knots), "Exposure_State")
  
  cv_fit <- cv.glmnet(design_matrix, data[["Mass"]], alpha = 1)
  
  m1 <- glmnet(design_matrix, data[["Mass"]], 
               alpha = 1, 
               lambda = 0.01)
  
  
  summary(m1)
  coef(m1)
  
  
  result = anova(model, model_reduced)
  aic[2] = AIC(model)
  loglik[2] = logLik(model)
  Test_statistic[2] = result$`F`[2]
  p_value[2] = result$`Pr(>F)`[2]
  
  data.table::data.table(Test = "S1",
                         State_1 = States[1],
                         State_2 = States[2],
                         Test_statistic = Test_statistic,
                         P_value = p_value,
                         Significant_difference = (p_value <= significance_level),
                         Time = NA,
                         Transformation = Transformation,
                         AIC = aic,
                         logLik = loglik)
}


### lasso

S_lasso <- function(data, significance_level = 0.05) {
  
}



### slope



S_slope <- function(data, significance_level = 0.05) {
  
}




