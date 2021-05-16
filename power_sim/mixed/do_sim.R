
source("./models.R")


library(powerHDX)
library(dplyr)
library(doParallel)
library(mgcv)
library(data.table)
library()


new_pf <- c(10, 15, 20, 30, 40, 50, 90, 100, 200, 210)
times = c(5, 10, 20, 30, 40, 50, 60, 100, 300, 500, 900, 1200, 1500, 1800,
          2100, 2400, 3600, 7200, 21600, 43200)

set.seed(17)
params <- readRDS("../all_params.RDS")
sequences <- sample(unique(params$sequence), 50, replace = FALSE)

all_params <- params %>% 
  dplyr::select(sequence, pH, step) %>% 
  unique() %>% 
  filter(sequence %in% sequences, pH == "6") %>% 
  group_by(sequence) %>% 
  slice(rep(1:n(), each = length(new_pf))) %>% 
  mutate(protection_factor = new_pf) %>% 
  ungroup() %>% 
  data.frame()


theoretical_spectra <- do.call(rbind, lapply(1:nrow(all_params), function(ith_row) {
  simulate_theoretical_spectra(sequence = all_params[ith_row, "sequence"],
                               charge = 1:3,
                               protection_factor = all_params[ith_row, "protection_factor"],
                               times = times,
                               pH = all_params[ith_row, "pH"],
                               temperature = 15,
                               n_molecules = 500,
                               time_step_const = 1,
                               use_markov = TRUE)
}))

spectra_by_seq <- split(theoretical_spectra, 
                        f = theoretical_spectra[, c("Sequence")])

print("Spectra ok!")

get_power = function(spectra_list) {
  lapply(
    spectra_list, function(spectrum) {
      tests_results = tryCatch({
        noisy_curves = get_noisy_deuteration_curves(spectrum, compare_pairs = TRUE,
                                                    reference = "all",
                                                    n_runs = 4,
                                                    n_replicates = 1,
                                                    mass_deviations = 5,
                                                    per_run_deviations = 0.1)
        calculate_hdx_power(noisy_curves,
                            tests = list(deuteros, Mixed_multi_knots_rep, Mixed_multi_knots_id),
                            significance_level  = 0.05, 
                            summarized = FALSE)
      }, error = function(e) {
        print(e)
        data.table::data.table()
        })
      seq = as.character(unique(noisy_curves[[1]][[1]]$Sequence))
      saveRDS(tests_results, file = paste("./results/", seq, "power.RDS", sep = "_"))
    }
  )
}

cores = detectCores()

get_power(spectra_by_seq)





