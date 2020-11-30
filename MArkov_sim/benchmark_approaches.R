library(dplyr)
library(powerHDX)
library(ggplot2)
library(microbenchmark)
library(profvis)
library(tidyverse)

source("working_one.R")


simulate_mean_exchange_for_aa = function(sequence, charge = NULL, protection_factor = 1,
                                         times = c(60, 600), pH = 7.5,
                                         temperature = 15, n_molecules = 100,
                                         time_step_const = 1, if_corr = 0,
                                         min_probability = 1e-4) {
  
  sequence = strsplit(sequence, "")[[1]]
  if (length(protection_factor) == 1L) {
    protection_factor = rep(protection_factor, length(sequence))
  }
  if (is.null(charge)) {
    charge = sample(2:6, 1)
  }
  
  peptide_iso_dist = get_approx_isotopic_distribution(sequence, min_probability)
  peptide_mass = peptide_iso_dist[[1]]
  isotopic_probs = peptide_iso_dist[[2]]
  maxND = peptide_iso_dist[[3]]
  maxD = peptide_iso_dist[[4]]
  kcHD = get_exchange_rates(sequence, "HD", pH, temperature, 'poly', if_corr)
  kcDH = get_exchange_rates(sequence, "DH", pH, temperature, 'poly')
  
  kmax = max(max(kcDH), max(kcHD))
  deltaT = time_step_const / kmax
  time_sequence = seq(0, max(times), deltaT)
  
  if (time_sequence == 0 && length(time_sequence) == 1) {
    print("There is no deuteration before given time point.")
    isotope_dists = data.frame()
  } else {
    
    transition_probs = get_exchange_probabilities(kcHD, kcDH, deltaT, protection_factor)
    transition_probs[["HH"]] <- 1 - transition_probs[["HD"]]
    transition_probs[["DD"]] <- 1 - transition_probs[["DH"]]
    
    HD_matrices_markov = get_HD_matrices_using_markov(sequence, transition_probs,
                                                      time_sequence, times,
                                                      n_molecules)
    
    HD_matrices_pkg = get_HD_matrices(sequence, transition_probs,
                                      time_sequence, times,
                                      n_molecules)
    
    exchange_pkg <- lapply(HD_matrices_pkg, colMeans)
    exchange_markov <- lapply(HD_matrices_markov, colMeans)
    
    
    result_pkg <- data.frame(Exchange = as.vector(do.call(c, exchange_pkg)),
                             AminoAcid = sequence,
                             Sequence = paste(sequence, collapse = ''),
                             Timepoint = sort(rep(times, length(sequence))),
                             Sim = "pkg")
    
    result_markov <- data.frame(Exchange = as.vector(do.call(c, exchange_markov)),
                                Sequence = paste(sequence, collapse = ''),
                                AminoAcid = sequence,
                                Timepoint = sort(rep(times, length(sequence))),
                                Sim = "markov")
    
    result <- rbind(result_pkg, result_markov)
    result$index <- rep(rep(1:length(sequence), length(times)), 2)
    result
  }
}


generate_curves_based_on_two_approaches <- function(one_param, n = 100) {
  result <- lapply(1:n, function(i) {
    rbind(curves_using_package(sequence = one_param$sequence,
                               protection_factor = one_param$protection_factor,
                               charge = one_param$charge,
                               times = times,
                               pH = one_param$pH,
                               time_step_const = one_param$step),
          curves_using_markov_chain(sequence = one_param$sequence,
                                    protection_factor = one_param$protection_factor,
                                    charge = one_param$charge,
                                    times = times,
                                    pH = one_param$pH,
                                    time_step_const = one_param$step))
  })
  do.call(rbind, result)
}

###### simulations


params <- readRDS(file = "../theo_spectra_sim/all_params_new_pf.RDS")
params$charge <- as.numeric(params$charge)
params$n_steps <- as.numeric(params$n_steps)
params$pH <- as.numeric(params$pH)
params$protection_factor <- as.numeric(params$protection_factor)


times <- c(5, 10, 20, 30, 40, 50, 60, 100, 300, 500, 900, 1200, 1500, 1800,
           2100, 2400, 3600, 7200, 21600, 43200)

set.seed(11)

example_params <- params %>% 
  filter(n_steps < 54805.0) %>% 
  sample_n(4, replace = FALSE)



### simulate matrices and compute mean exchange for every amino acid from the sequence

example_param = example_params[3, ]

split_sequence <- strsplit(example_param$sequence, "")[[1]]

result_exchange <- simulate_mean_exchange_for_aa(sequence = example_param$sequence,
                                                 protection_factor = example_param$protection_factor,
                                                 charge = example_param$charge,
                                                 times = times,
                                                 pH = example_param$pH,
                                                 time_step_const = example_param$step)


mean_exchange <- ggplot(result_exchange, aes(x = as.character(index), y = Exchange , fill = Sim)) +
  geom_bar(stat = 'identity', position = "dodge") +
  scale_x_discrete(labels = split_sequence) +
  facet_wrap(~Timepoint) +
  ggtitle(paste(split_sequence, collapse = ''))

mean_exchange_point_line <- ggplot(result_exchange, aes(x = index, 
                                                        y = Exchange, 
                                                        col = Sim)) +
  geom_point() +
  geom_line() +
  scale_x_discrete(labels = split_sequence) +
  facet_wrap(~Timepoint) +
  ggtitle(paste(split_sequence, collapse = ''))



ggsave("mean_exchange_from_matrices.png", 
       plot = mean_exchange,
       device = "png")

ggsave("mean_exchange_point_line.png", 
       plot = mean_exchange_point_line,
       device = "png")



### generate curves

generated_curves_100 <- do.call(rbind, lapply(1:nrow(example_params), function(i) {
  generate_curves_based_on_two_approaches(example_params[i, ], n = 100)
}))

saveRDS(generated_curves_100, file = "./generated_curves_100.RDS")

result <- readRDS("G:/HDX_sim_vol2/MArkov_sim/generated_curves_100.RDS")


#####

ggplot(result, aes(x = as.factor(Exposure), y = Mass, col = Sim)) +
  geom_boxplot(position = "dodge") +
  facet_wrap(~Sequence, scales = "free_y") 


ggplot(result[Sequence == "FPTTKTY", ], aes(x = as.factor(Exposure), y = Mass, col = Sim)) +
  geom_boxplot(position = "dodge") +
  geom_jitter(shape=16, position=position_jitter(0.2))


#### average

average_curves <- result %>% 
  group_by(Exposure, Sim, Sequence, PF, PH) %>% 
  mutate(Mean_mass = mean(Mass)) %>% 
  ungroup %>% 
  select(-Mass) %>% 
  unique() %>% 
  rename(Mass = Mean_mass) %>% 
  group_by(Exposure, Sim, Sequence, PF, PH)


ggplot(average_curves, aes(x = Exposure, y = Mass, col = Sim)) +
  geom_point() +
  geom_line() +
  facet_wrap(~Sequence, scales = "free_y")


#######


diff <- average_curves %>% 
  group_by(Exposure, Sequence) %>% 
  summarise(diff = abs(Mass[Sim == "pkg"] - Mass[Sim == "markov"]))

ggplot(diff, aes(x = Exposure, y = diff, col = Sequence)) +
  geom_line() +
  geom_point() +
  ggtitle("|Mass(sim) - Mass(markov)|") +
  geom_hline(yintercept = 0, size = 1)


