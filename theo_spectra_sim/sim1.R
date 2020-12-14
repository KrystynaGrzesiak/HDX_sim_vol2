library(doParallel)
library(powerHDX)

setwd("~/HDX_sim_vol2/theo_spectra_sim")
load("sysdata.rda")

sim_theo_spectra = function(all_params, times) {
  lapply(1:nrow(all_params), function(ith_row) {
    print(paste("Simulation", ith_row, "\n"))
    res = tryCatch(simulate_theoretical_spectra(sequence = all_params[ith_row, "sequence"],
                                                charge = all_params[ith_row, "charge"],
                                                protection_factor = all_params[ith_row, "protection_factor"],
                                                times = times,
                                                pH = all_params[ith_row, "pH"],
                                                temperature = 15,
                                                n_molecules = 500,
                                                time_step_const = all_params[ith_row, "step"],
                                                use_markov = TRUE),
                   error = function(e) {
                     print(e)
                     data.frame()})
    saveRDS(res, file = paste0("./results/theo_spectrum_", all_params[ith_row, "sequence"], "_",
                               all_params[ith_row, "pH"], "_",
                               all_params[ith_row, "protection_factor"], "_",
                               ith_row, ".RDS", collapse = ""))
    rm(res)
  })
}

times = c(5, 10, 20, 30, 40, 50, 60, 100, 300, 500, 900, 1200, 1500, 1800,
          2100, 2400, 3600, 7200, 21600, 43200)


all_params = readRDS("./all_params_new_pf.RDS")
all_params$sequence = as.character(all_params$sequence)
all_params$charge = as.numeric(as.character(all_params$charge))
set.seed(1410)

sim_theo_spectra(all_params, times)



