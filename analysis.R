library(powerHDX)
library(dplyr)
library(ggplot2)

# reading spectra -----------------------------

files <- list.files(path = "./theo_spectra_sim/results", 
                    pattern = "RDS$", full.names = TRUE)

sim_results <- do.call("rbind", lapply(files, readRDS))

spectra_by_ph_sequence <- split(sim_results, f = sim_results[, c("Sequence", "PH")])

saveRDS(spectra_by_ph_sequence, file = "spectra_by_ph_seq.RDS")

theo_deut_curves <- lapply(spectra_by_ph_seq, function(spectrum) {
  get_deuteration_curve_single_spectrum(spectrum)
})

saveRDS(theo_deut_curves, file = "theo_deut_curves.RDS")

crv <- do.call(rbind, readRDS("./theo_spectra_sim/theo_deut_curves.RDS"))


p <- ggplot(crv, aes(x = Exposure, y = Mass, color = factor(PF),
                shape = factor(PH), linetype = factor(PH))) +
  geom_point() +
  geom_line() +
  scale_color_discrete("PF") +
  theme_bw() +
  theme(legend.position = "bottom") +
  facet_wrap(~ Sequence, scales = "free_y")

cairo_pdf("./results/curves-without-noise.pdf", width = 50, height = 40)
p
dev.off()

# adding noise -----------------

