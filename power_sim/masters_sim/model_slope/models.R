library(lmerTest)
library(glmnet)

truncated_lines <- function(x, knots){
  sapply(knots, function(kappa) {
    (x - kappa)*(x > kappa)
  })
}


houde <- function(data, significance_level = 0.05) {
  
  States = unique(data$State)
  confidence_limit = 1 - significance_level
  
  t_value <- qt(c(significance_level/2, 1 - significance_level/2), df = 2)[2]
  
  calc_dat <- data %>%
    group_by(Sequence, State, Exposure, Rep, Experimental_state) %>%
    summarise(avg_exp_mass = mean(Mass)) %>% #usredniamy po ladunkach
    ungroup() %>%
    group_by(State, Sequence, Exposure, Experimental_state) %>%
    summarize(deut_uptake = mean(avg_exp_mass),
              err_avg_mass = sd(avg_exp_mass)/sqrt(length(Exposure))) %>% #usredniamy po replikacjach i liczymy niepewnosc
    group_by(State, Sequence, Experimental_state) %>%
    mutate(err_deut_uptake = sqrt(err_avg_mass^2 + err_avg_mass[Exposure == 0]^2)) %>%
    ungroup(.) %>%
    group_by(Sequence, Exposure) %>%
    summarise(diff_deut_uptake = deut_uptake[Experimental_state == "A"] -
                deut_uptake[Experimental_state == "B"],
              err_diff_deut_uptake = sqrt(err_deut_uptake[Experimental_state == "A"]^2 +
                                            err_deut_uptake[Experimental_state == "B"]^2)) #liczymy roznice w deuterium uptake miedzy stanami
  
  
  avg_difference <- mean(calc_dat$diff_deut_uptake)
  
  x_threshold <- t_value * mean(calc_dat[["err_diff_deut_uptake"]], na.rm = TRUE)/sqrt(length(calc_dat))
  
  data.table(Test = "Houde",
             State_1 = States[1],
             State_2 = States[2],
             Test_statistic = NA,
             P_value = NA,
             Significant_difference = abs(avg_difference) > x_threshold,
             Time = NA,
             Transformation = NA,
             AIC = NA,
             logLik = NA)
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
    model <- lmer(Mass ~ Exposure*Experimental_state + (1|Exposure) + (Exposure || Experimental_state) + X_reduced,
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
  
  
  data.table(Test = "S_model_random_state_exposure",
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