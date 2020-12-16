library(powerHDX)
library(ggplot2)

load("sysdata.rda")

times = c(5, 10, 20, 30, 40, 50, 60, 100, 300, 500, 900, 1200, 1500, 1800,
          2100, 2400, 3600, 7200, 21600, 43200)


all_params = readRDS("./all_params_new_pf.RDS")
all_params$sequence = as.character(all_params$sequence)
all_params$charge = as.numeric(as.character(all_params$charge))
set.seed(1410)

# spec <- readRDS("G:/HDX_sim_vol2/theo_spectra_sim/results/theo_spectrum_HHFGKEFTPPVQAA_7.5_100_8035.RDS")

huge_vec_param <- all_params[8035,]

huge_vec_param <- all_params[8040,]


simulate_theoretical_spectra(sequence = huge_vec_param[["sequence"]],
                             charge = huge_vec_param[["charge"]],
                             protection_factor = huge_vec_param[["protection_factor"]],
                             times = times,
                             pH = huge_vec_param[["pH"]],
                             temperature = 15,
                             n_molecules = 500,
                             time_step_const = huge_vec_param[["step"]],
                             use_markov = TRUE)





