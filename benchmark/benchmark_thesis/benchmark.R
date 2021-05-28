library(powerHDX)
library(microbenchmark)
library(dplyr)

simulate_markov <- function(params) {
  simulate_theoretical_spectra(sequence = params[["sequence"]],
                               charge = 3,
                               protection_factor = 100,
                               times = 43200,
                               pH = 7,
                               temperature = 15,
                               n_molecules = 500,
                               time_step_const = 1,
                               use_markov = TRUE)
}

simulate_rcpp <- function(params) {
  simulate_theoretical_spectra(sequence = params[["sequence"]],
                               charge = 3,
                               protection_factor = 100,
                               times = 43200,
                               pH = 7,
                               temperature = 15,
                               n_molecules = 500,
                               time_step_const = 1,
                               use_markov = FALSE)
}

new_pf <- c(100)
set.seed(10)
params <- readRDS("./all_params.RDS")

sample_data <- params %>% 
  select(sequence) %>% 
  unique() %>% 
  mutate(lengths = nchar(sequence)) %>% 
  mutate(intervals = cut(lengths, seq(0, 60, 15))) %>% 
  group_by(intervals) %>% 
  sample_n(3) %>% 
  ungroup() %>% 
  arrange(lengths)


lapply(1:nrow(sample_data), function(i) {
  benchmark_markov <- microbenchmark({simulate_markov(sample_data[i, ])},
                                     times = 10)
  benchmark_rcpp <- microbenchmark({simulate_rcpp(sample_data[i, ])},
                                   times = 10)
  
  res <- data.frame(markov_times = benchmark_markov$time,
                    rcpp_times = benchmark_rcpp$time,
                    sample_data[i, ])
  
  saveRDS(res, paste0("./results/benchmark_", sample_data[i, "sequence"], ".RDS"))
})


