
library(powerHDX)
library(dplyr)
library(ggplot2)

files <- list.files(path = "G:/HDX_sim_vol2/theo_spectra_sim/results", pattern = "\\.RDS$", full.names = TRUE) 
sim_results <- do.call("rbind", lapply(files, readRDS))

params_from_file <- readRDS("all_params_new_pf.RDS") 

nrow(params_from_file) #all spectra #12117 rows
length(files) #spectra from simulation #11117 files
nrow(params_from_file) - length(files) #missing spectra

missing1 <- sum(sapply(files, function(i) {
  dat = readRDS(i)
  nrow(dat) == 0
}))

missing1 #empty data frame, i.e. error during simulation #18 data frames


params <- params_from_file %>% 
  select(sequence, pH, charge, protection_factor) %>% 
  unique() %>% 
  nrow()

colnames(params) <- c("Sequence", "PH", "Charge", "PF")
params[["Charge"]] <- as.numeric(params[["Charge"]])



sim_params <- sim_results %>% 
  select(Sequence, PH, Charge, PF) %>% 
  unique() 

missing_params <- setdiff(params, sim_params) 
colnames(missing_params) <- c("sequence", "pH", "charge", "protection_factor")

saveRDS(missing_params, file = "missing_params.RDS")


# spectra_by_ph_sequence = split(sim_results, f = sim_results[, c("Sequence", "PH")])
# saveRDS(spectra_by_ph_sequence, file = "spectra_by_ph_seq.RDS")
# 
# theo_deut_curves = lapply(spectra_by_ph_seq, function(spectrum) {
#   get_deuteration_curve_single_spectrum(spectrum)
# })
# 
# saveRDS(theo_deut_curves, file = "theo_deut_curves.RDS")
