
library(doParallel)
library(powerHDX)

setwd("~/HDX_sim_vol2/theo_spectra_sim")
load("sysdata.rda")

sim_theo_spectra = function(all_params, n_cpus, times) {
  all_params$cpu = sample(1:n_cpus, nrow(all_params), replace = TRUE)
  all_params = split(all_params, all_params$cpu)
  
  mclapply(all_params, function(df_per_cpu) {
    df_per_cpu$sequence = as.character(df_per_cpu$sequence)
    df_per_cpu$charge = as.numeric(as.character(df_per_cpu$charge))
    lapply(1:nrow(df_per_cpu), function(ith_row) {
      print(paste("Simulation", ith_row, "\n"))
      res = tryCatch(simulate_theoretical_spectra(sequence = df_per_cpu[ith_row, "sequence"],
                                                  charge = df_per_cpu[ith_row, "charge"],
                                                  protection_factor = df_per_cpu[ith_row, "protection_factor"],
                                                  times = times,
                                                  pH = df_per_cpu[ith_row, "pH"],
                                                  temperature = 15,
                                                  n_molecules = 500,
                                                  time_step_const = df_per_cpu[ith_row, "step"],
                                                  use_markov = TRUE),
                     error = function(e) {
                       print(e)
                       data.frame()})
      saveRDS(res, file = paste0("./missing_results/theo_spectrum_", df_per_cpu[ith_row, "sequence"], "_",
                                 df_per_cpu[ith_row, "pH"], "_",
                                 df_per_cpu[ith_row, "protection_factor"], "_",
                                 ith_row, ".RDS", collapse = ""))
      rm(res)
    })
  }, mc.cores = n_cpus)
}

times = c(5, 10, 20, 30, 40, 50, 60, 100, 300, 500, 900, 1200, 1500, 1800,
          2100, 2400, 3600, 7200, 21600, 43200)


all_params = readRDS("../theo_spectra_sim/missing_params.RDS")
all_params$sequence = as.character(all_params$sequence)
all_params$charge = as.numeric(as.character(all_params$charge))

cores = detectCores()

sim_theo_spectra(all_params, cores, times)













