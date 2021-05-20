

truncated_lines <- function(x, kappa){
  (x - kappa)*(x > kappa)
}

Mixed_double_knots_id_3600 <- function(data, significance_level = 0.05) {
  
  data[["id"]] <- paste0(data$Rep, data$Charge, data$Experimental_state)
  
  transformed_data <- data
  States = unique(data$State)
  
  knots <- unique(setdiff(transformed_data$Exposure, c(3600, max(transformed_data$Exposure), min(transformed_data$Exposure))))
  Test = aic = loglik = Test_statistic = p_value = rep(NA, length(knots))
  
  data.table::data.table(do.call(rbind, lapply(1:length(knots), function(knot) {
    
    knot_points <- unique(c(knots[knot], 3600))
    
    X <- sapply(1:length(knot_points), function(k) {
      truncated_lines(transformed_data$Exposure, knot_points[k])
    })
    
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
    
    
    Test = paste0(knots[knot], "_and_3600", "_id")
    
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
  })))
}
