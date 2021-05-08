

truncated_lines <- function(x, kappa){
  (x - kappa)*(x > kappa)
}



S_multi_knots <- function(data, significance_level = 0.05) {
  States = unique(data$State)
  Transformation = c("identity", "log")
  data[["Exposure"]] = log(data[["Exposure"]] + 1)
  
  knots <- unique(setdiff(data$Exposure, c(max(data$Exposure), min(data$Exposure))))
  Test = aic = loglik = Test_statistic = p_value = rep(NA, length(knots))

  data.table::data.table(do.call(rbind, lapply(1:length(knots), function(knot) {
    
    X <- sapply(1:length(knot), function(k) {
      truncated_lines(data$Exposure, knot[k])
    })
    
    model = lm(Mass ~ State*Exposure + X,
               data = data)
    model_reduced = lm(Mass ~ Exposure + X,
                       data = data)
    result = anova(model, model_reduced)
    aic = AIC(model)
    loglik = logLik(model)
    Test_statistic = result$`F`[2]
    p_value = result$`Pr(>F)`[2]
    Test = paste0("Spline_", exp(knots[knot])-1)
    
    data.table::data.table(Test = Test,
                           State_1 = States[1],
                           State_2 = States[2],
                           Test_statistic = Test_statistic,
                           P_value = p_value,
                           Significant_difference = (p_value <= significance_level),
                           Time = NA,
                           Transformation = "log",
                           AIC = aic,
                           logLik = loglik)
  })))
}



