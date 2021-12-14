library(powerHaDeX)

times = c(0.06, 10, 60, 300, 1500, 7200, 86400)

spec1 <- simulate_theoretical_spectra("PPAQHI", protection_factor = 10, charge = 1:3, times = times)
spec2 <- simulate_theoretical_spectra("PPAQHI", protection_factor = 20, charge = 1:3, times = times)
spec3 <- simulate_theoretical_spectra("PPAQHI", protection_factor = 30, charge = 1:3, times = times)

spectra <- rbind(spec1, spec2, spec3)

a <- get_noisy_deuteration_curves(spectra, 
                                  compare_pairs = TRUE, 
                                  reference = "all", 
                                  n_replicates = 4,
                                  n_experiments  = 10)
data <- a[[6]][[10]]




a = calculate_hdx_power(data, tests = list(S_lasso_ridge, deuteros), summarized = FALSE)

a %>% 
  filter(State_1 ==10, State_2 == 20, Transformation == "identity")
