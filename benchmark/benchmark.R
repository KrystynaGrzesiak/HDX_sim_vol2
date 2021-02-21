
library(powerHDX)
library(microbenchmark)
library(dplyr)

new_pf <- c(100)
set.seed(10)
params <- readRDS("./all_params.RDS")
sequences <- sample(unique(params$sequence), 100, replace = FALSE)

all_params <- params %>% 
  dplyr::select(sequence, pH, step) %>% 
  unique() %>% 
  filter(sequence %in% sequences, pH == "6") %>% 
  group_by(sequence) %>% 
  slice(rep(1:n(), each = length(new_pf))) %>% 
  mutate(protection_factor = new_pf) %>% 
  mutate(charge = sample(1:5,length(new_pf),  replace = TRUE)) %>% 
  ungroup() %>% 
  data.frame()


simulate_markov <- function(params, times) {
  simulate_theoretical_spectra(sequence = params[["sequence"]],
                               charge = params[["charge"]],
                               protection_factor = params[["protection_factor"]],
                               times = times,
                               pH = params[["pH"]],
                               temperature = 15,
                               n_molecules = 500,
                               time_step_const = params[["step"]],
                               use_markov = TRUE)
}
  
simulate_rcpp <- function(params, times) {
  simulate_theoretical_spectra(sequence = params[["sequence"]],
                               charge = params[["charge"]],
                               protection_factor = params[["protection_factor"]],
                               times = times,
                               pH = params[["pH"]],
                               temperature = 15,
                               n_molecules = 500,
                               time_step_const = params[["step"]],
                               use_markov = FALSE)
}
  
times = c(5, 10, 20, 30, 40, 50, 60, 100, 300, 500, 900, 1200, 1500, 1800,
          2100, 2400, 3600, 7200, 21600, 43200)

lapply(1:nrow(all_params), function(i) {
  benchmark <- microbenchmark({simulate_rcpp(all_params[i, ], times)}, 
                              {simulate_markov(all_params[i, ], times)},
                              times = 10)
  saveRDS(benchmark, paste0("./results/benchmark_", i, ".RDS"))
})





