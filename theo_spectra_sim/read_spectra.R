
library(powerHDX)
library(dplyr)
library(ggplot2)

files <- list.files(path = "G:/HDX_sim_vol2/theo_spectra_sim/results", pattern = "\\.RDS$", full.names = TRUE)
sim_results <- do.call("rbind", lapply(files, readRDS))

spectra_by_ph_sequence = split(sim_results, f = sim_results[, c("Sequence", "PH")])
saveRDS(spectra_by_ph_sequence, file = "spectra_by_ph_seq.RDS")

theo_deut_curves = lapply(spectra_by_ph_seq, function(spectrum) {
  get_deuteration_curve_single_spectrum(spectrum)
})

saveRDS(theo_deut_curves, file = "theo_deut_curves.RDS")
