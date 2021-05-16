

truncated_lines <- function(x, kappa){
  (x - kappa)*(x > kappa)
}



S_multi_knots <- function(data, significance_level = 0.05) {
  
  Transformation = c("log", "identity")
  data.table::rbindlist(lapply(Transformation, function(trans) {
    
    if(trans == "log") {
      transformed_data <- data
      transformed_data[["Exposure"]] = log(transformed_data[["Exposure"]] + 1)
    } else {
      transformed_data <- data
    }
    
    States = unique(data$State)
    
    knots <- unique(setdiff(transformed_data$Exposure, c(max(transformed_data$Exposure), min(transformed_data$Exposure))))
    Test = aic = loglik = Test_statistic = p_value = rep(NA, length(knots))
    
    data.table::data.table(do.call(rbind, lapply(1:length(knots), function(knot) {
      X <- sapply(1:length(knot), function(k) {
        truncated_lines(transformed_data$Exposure, knot[k])
      })
      
      model = lm(Mass ~ State*Exposure + X,
                 data = transformed_data)
      model_reduced = lm(Mass ~ Exposure + X,
                         data = transformed_data)
      result = anova(model, model_reduced)
      aic = AIC(model)
      loglik = as.numeric(logLik(model))
      Test_statistic = result$`F`[2]
      p_value = result$`Pr(>F)`[2]
      
      if(trans == "log") {
        Test = paste0("Spline_", exp(knots[knot])-1)
      } else {
        Test = knots[knot]
      }
      
      data.table::data.table(Test = Test,
                             State_1 = States[1],
                             State_2 = States[2],
                             Test_statistic = Test_statistic,
                             P_value = p_value,
                             Significant_difference = (p_value <= significance_level),
                             Time = NA,
                             Transformation = trans,
                             AIC = aic,
                             logLik = loglik)
    })))
  }))
}



