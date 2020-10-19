
library(powerHDX)
library(dplyr)
library(ggplot2)
library(plyr)

files <- list.files(path = "G:/HDX_sim_vol2/theo_spectra_sim/results", pattern = "\\.RDS$", full.names = TRUE) 
sim_results <- do.call("rbind", lapply(files, readRDS))

params_from_file <- readRDS("all_params_new_pf.RDS") 

nrow(params_from_file) #all spectra #12 117 rows
length(files) #spectra from simulation #11 107 files
nrow(params_from_file) - length(files) #missing spectra

missing1 <- sum(sapply(files, function(i) {
  dat = readRDS(i)
  nrow(dat) == 0
}))

missing1 #empty data frame, i.e. error during simulation #18 data frames


params_from_file[["charge"]] <- as.numeric(params_from_file[["charge"]]) 

params <- params_from_file %>% 
  select(-step, -n_steps, -time)


sim_params <- sim_results %>% 
  select(Sequence, PH, Charge, PF)

colnames(sim_params) <- c("sequence", "pH", "charge", "protection_factor")


sim_params_step <- merge(x = sim_params, y = params, 
                         by = c("sequence", "pH", "charge", "protection_factor"), 
                         all.x = TRUE) %>% 
  unique()



missing_params <- setdiff(params, sim_params_step) #4131 rows
colnames(missing_params) <- c("sequence", "pH", "charge", "protection_factor", "size_of_time_step")

saveRDS(missing_params, file = "missing_params.RDS")


# spectra_by_ph_sequence = split(sim_results, f = sim_results[, c("Sequence", "PH")])
# saveRDS(spectra_by_ph_sequence, file = "spectra_by_ph_seq.RDS")
# 
# theo_deut_curves = lapply(spectra_by_ph_seq, function(spectrum) {
#   get_deuteration_curve_single_spectrum(spectrum)
# })
# 
# saveRDS(theo_deut_curves, file = "theo_deut_curves.RDS")
