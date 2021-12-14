
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

