
library(powerHDX)
library(dplyr)
library(ggplot2)
library(plyr)

params_from_file <- readRDS("all_params_new_pf.RDS") 

ith_row <- 1

times <- c(5, 10, 20, 30, 40, 50, 60, 100, 300, 500, 900, 1200, 1500, 1800,
           2100, 2400, 3600, 7200, 21600, 43200)

param_big_vec <- params_from_file %>% 
  filter(sequence == "SALSDLHAHKLRVDPV", 
         pH == 7.5, protection_factor == 10, 
         charge == 4)


param_big_vec$sequence <- as.character(param_big_vec$sequence)
param_big_vec$charge <- as.numeric(as.character(param_big_vec$charge))

simulate_theoretical_spectra(sequence = param_big_vec[ith_row, "sequence"],
                             charge = param_big_vec[ith_row, "charge"],
                             protection_factor = param_big_vec[ith_row, "protection_factor"],
                             times = times,
                             pH = param_big_vec[ith_row, "pH"],
                             temperature = 15,
                             n_molecules = 500,
                             time_step_const = param_big_vec[ith_row, "step"])


simulate_theoretical_spectra(sequence = params_from_file[ith_row, "sequence"],
                             charge = params_from_file[ith_row, "charge"],
                             protection_factor = params_from_file[ith_row, "protection_factor"],
                             times = times,
                             pH = params_from_file[ith_row, "pH"],
                             temperature = 15,
                             n_molecules = 500,
                             time_step_const = params_from_file[ith_row, "step"])




x <- runif(100000000)

y <- seq(0.1, 1, length.out = 20)

separated_times <- split(x, cut(x, c(0, y), include.lowest = TRUE))


data.frame(time = x, a = 0, b = x)











