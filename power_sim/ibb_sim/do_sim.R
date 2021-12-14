source("./models.R")


library(powerHaDeX)
library(dplyr)
library(doParallel)
library(mgcv)
library(data.table)

tests_fun[[4]] <- test_semiparametric
tests_fun[[5]] <- test_houde


new_pf <- c(10, 15, 20, 30, 40, 50, 90, 100, 110, 190, 200, 210, 300, 400)


times = c(0.06, 10, 60, 300, 1500, 7200, 86400)

set.seed(17)
params <- readRDS("../all_params.RDS")
sequences <- sample(unique(params$sequence), 100, replace = FALSE)

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
        browser()
        noisy_curves = get_noisy_deuteration_curves(spectrum, compare_pairs = TRUE,
                                                    reference = "all",
                                                    n_replicates = 4,
                                                    n_experiments = 100,
                                                    per_run_deviations = runif(21, 0, 0.5), 
                                                    mass_deviations = rnorm(21, 50, 10))
        calculate_hdx_power(noisy_curves,
                            tests = tests_fun,
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
