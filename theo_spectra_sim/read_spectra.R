
library(powerHDX)
library(dplyr)
library(ggplot2)

files <- list.files(path = "G:/HDX_sim_vol2/theo_spectra_sim/results", pattern = "\\.RDS$", full.names = TRUE)
sim_results <- do.call("rbind", lapply(files, readRDS))

params <- readRDS("all_params_new_pf.RDS") %>% 
  select(sequence, pH, protection_factor, charge) %>% 
  mutate(Charge = as.numeric(charge)) %>% 
  select(-charge) %>% 
  unique()

colnames(params) = c("Sequence", "PH", "PF", "Charge")

sim_params <- sim_results %>% 
  select(Sequence, PH, PF, Charge) %>% 
  unique()

setdiff(params, sim_params) %>% 
  nrow()




spectra_by_ph_sequence = split(sim_results, f = sim_results[, c("Sequence", "PH")])
saveRDS(spectra_by_ph_sequence, file = "spectra_by_ph_seq.RDS")

theo_deut_curves = lapply(spectra_by_ph_seq, function(spectrum) {
  get_deuteration_curve_single_spectrum(spectrum)
})

saveRDS(theo_deut_curves, file = "theo_deut_curves.RDS")
