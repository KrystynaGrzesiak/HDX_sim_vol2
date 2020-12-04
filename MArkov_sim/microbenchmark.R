library(powerHDX)
library(ggplot2)
library(dplyr)
library(microbenchmark)

load("../theo_spectra_sim/sysdata.rda")


###### simulations


params <- readRDS(file = "../theo_spectra_sim/all_params_new_pf.RDS")
params$charge <- as.numeric(params$charge)
params$n_steps <- as.numeric(params$n_steps)
params$pH <- as.numeric(params$pH)
params$protection_factor <- as.numeric(params$protection_factor)


times <- c(5, 10, 20, 30, 40, 50, 60, 100, 300, 500, 900, 1200, 1500, 1800,
           2100, 2400, 3600, 7200, 21600, 43200)

set.seed(17)

example_params <- params %>% 
  sample_n(1, replace = FALSE)




benchmark = microbenchmark({
  lapply(1:nrow(example_params), function(i) {
    simulate_theoretical_spectra(sequence = example_params$sequence[i],
                                 protection_factor = example_params$protection_factor[i],
                                 charge = example_params$charge[i],
                                 times = times,
                                 pH = example_params$pH[i],
                                 time_step_const = example_params$step[i],
                                 use_markov = TRUE)
  })
}, {
  lapply(1:nrow(example_params), function(i) {
    simulate_theoretical_spectra(sequence = example_params$sequence[i],
                                 protection_factor = example_params$protection_factor[i],
                                 charge = example_params$charge[i],
                                 times = times,
                                 pH = example_params$pH[i],
                                 time_step_const = example_params$step[i],
                                 use_markov = FALSE)
  })
}, times = 20L)




saveRDS(benchmark, "microbench.RDS")



##########

example_params <- params %>% 
  filter(sequence == "TPAVH", pH == 6) %>% 
  sample_n(1, replace = FALSE)


benchmark2 = microbenchmark({
  lapply(1:nrow(example_params), function(i) {
    simulate_theoretical_spectra(sequence = example_params$sequence[i],
                                 protection_factor = example_params$protection_factor[i],
                                 charge = example_params$charge[i],
                                 times = times,
                                 pH = example_params$pH[i],
                                 time_step_const = example_params$step[i],
                                 use_markov = TRUE)
  })
}, {
  lapply(1:nrow(example_params), function(i) {
    simulate_theoretical_spectra(sequence = example_params$sequence[i],
                                 protection_factor = example_params$protection_factor[i],
                                 charge = example_params$charge[i],
                                 times = times,
                                 pH = example_params$pH[i],
                                 time_step_const = example_params$step[i],
                                 use_markov = FALSE)
  })
}, times = 20L)




saveRDS(benchmark2, "microbench1.RDS")










