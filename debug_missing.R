library(powerHDX)
library(dplyr)
library(ggplot2)
library(plyr)


times <- c(5, 10, 20, 30, 40, 50, 60, 100, 300, 500, 900, 1200, 1500, 1800,
          2100, 2400, 3600, 7200, 21600, 43200)

params_from_file <- readRDS("all_params_new_pf.RDS") 
missing_params <- readRDS("missing_params.RDS")



little_sim <- function(all_params, times) {
  all_params$sequence <- as.character(all_params$sequence)
  all_params$charge <- as.numeric(as.character(all_params$charge))
  lapply(1:nrow(all_params), function(ith_row) {
    simulate_theoretical_spectra(sequence = all_params[ith_row, "sequence"],
                                 charge = all_params[ith_row, "charge"],
                                 protection_factor = all_params[ith_row, "protection_factor"],
                                 times = times,
                                 pH = all_params[ith_row, "pH"],
                                 temperature = 15,
                                 n_molecules = 500,
                                 time_step_const = all_params[ith_row, "step"])
  })
}


little_sim(missing_params[1:100, ], c(5, 10, 20, 30, 40, 50, 60, 100, 300))


