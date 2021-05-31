

truncated_lines <- function(x, kappa){
  (x - kappa)*(x > kappa)
}



Mixed_multi_knots_rep <- function(data, significance_level = 0.05) {
  
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
      
      # continuous, identity
      model = lmerTest::lmer(Mass ~ Exposure*State + (1|Rep) + X,
                             data = transformed_data,
                             REML = FALSE)
      model_reduced = lmerTest::lmer(Mass ~ Exposure + (1|Rep) + X,
                                     data = transformed_data,
                                     REML = FALSE)
      result = anova(model, model_reduced)
      aic = AIC(model)
      loglik = as.numeric(logLik(model))
      Test_statistic = result$Chisq[2]
      p_value = result$`Pr(>Chisq)`[2]
      
      
      if(trans == "log") {
        Test = paste0("Spline_", exp(knots[knot])-1, "_rep")
      } else {
        Test = paste0(knots[knot], "_rep")
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




Mixed_multi_knots_id <- function(data, significance_level = 0.05) {
  
  Transformation = c("log", "identity")
  
  data[["id"]] <- paste0(data$Rep, data$Charge, data$Experimental_state)
  
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
      
      X <- truncated_lines(transformed_data$Exposure, knots[knot])
      
      # continuous, identity
      model = lmerTest::lmer(Mass ~ Exposure*State + (1|id) + X,
                             data = transformed_data,
                             REML = FALSE)
      model_reduced = lmerTest::lmer(Mass ~ Exposure + (1|id) + X,
                                     data = transformed_data,
                                     REML = FALSE)
      result = anova(model, model_reduced)
      aic = AIC(model)
      loglik = as.numeric(logLik(model))
      Test_statistic = result$Chisq[2]
      p_value = result$`Pr(>Chisq)`[2]
      
      
      if(trans == "log") {
        Test = paste0("Spline_", exp(knots[knot])-1, "_id")
      } else {
        Test = paste0(knots[knot], "_id")
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




