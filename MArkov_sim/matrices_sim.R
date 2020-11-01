library(powerHDX)
library(ggplot2)
library(dplyr)
library(markovchain)

load("../theo_spectra_sim/sysdata.rda")

simulate_states_dist = function(sequence, charge = NULL, protection_factor = 1,
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
    times_to_record = get_recording_times(time_sequence, times)
    times_to_record = setdiff(times_to_record, 0)
    
    transition_probs = get_exchange_probabilities(kcHD, kcDH, deltaT, protection_factor)
    HD_matrices = get_HD_matrices(sequence, transition_probs,
                                  time_sequence, times_to_record,
                                  n_molecules)
    list(matrices = HD_matrices, 
         probabilities = transition_probs, 
         time_sequences = split(time_sequence,
                                cut(time_sequence, c(0, times_to_record),
                                    include.lowest = TRUE)))
  }
}

#compute probabilities and simulate matrices

params <- readRDS(file = "../theo_spectra_sim/all_params_new_pf.RDS")
param <- params[117, ]
times <- c(5, 10, 20, 30, 40, 50, 60, 100, 300, 500, 900, 1200, 1500, 1800,
           2100, 2400, 3600, 7200, 21600, 43200)


set.seed(17)

matrices_and_probs <- simulate_states_dist(sequence = param$sequence,
                                           protection_factor = param$protection_factor,
                                           charge = param$charge,
                                           times = times,
                                           pH = param$pH,
                                           time_step_const = param$step)


### package simulation (HD matrices)

split_sequence <- strsplit(param$sequence, "")[[1]]

matrices_from_sim <- matrices_and_probs[["matrices"]]
names(matrices_from_sim) <- times

mean_exchange <- lapply(1:length(times), function(i) {
  df <- data.frame(t(rbind(colMeans(matrices_from_sim[[i]]), times[i], split_sequence)))
  names(df) <- c("mean_exchange", "time", "amino_acid")
  df$mean_exchange <- as.numeric(df$mean_exchange)
  df$time <- as.numeric(df$time)
  df$unique_acid <- paste0(1:length(df$amino_acid), df$amino_acid)
  df
})


result <- do.call(rbind, mean_exchange)


ggplot(result, aes(x = unique_acid, y = mean_exchange)) +
  geom_bar(stat = 'identity') + 
  scale_x_discrete(labels = split_sequence) +
  facet_wrap(~result$time)+
  ylim(0, 0.3)



### Markov chain approach

transition_probs <- matrices_and_probs[["probabilities"]]

transition_probs[["HH"]] <- 1 - transition_probs[["HD"]]
transition_probs[["DD"]] <- 1 - transition_probs[["DH"]]

time_sequences <- matrices_and_probs[["time_sequences"]]
names(time_sequences) <- times

#transition matrices for Markov chain


transition_matrices <- lapply(1:length(split_sequence), function(i) {
  matrix(c(transition_probs[["HH"]][i], transition_probs[["HD"]][i],
           transition_probs[["DH"]][i], transition_probs[["DD"]][i]), nrow = 2, 
         byrow = TRUE)
})


chain1 <- new('markovchain', 
              states = c("H", "D"), 
              transitionMatrix = transition_matrices[[1]],
              name = "Deuteration")

initial_state <- c(1, 0)



mean_exchange_markov <- lapply(1:length(times), function(time) {
  mean_exchange = sapply(1:length(split_sequence), function(amino_acid) {
    chain <- new('markovchain', states = c("0", "1"), 
                 transitionMatrix = transition_matrices[[amino_acid]])
    mean(sample(c(0, 1), 100, replace = TRUE, prob = initial_state*(chain^length(time_sequences[[time]]))))
  })
  data.frame(mean_exchange = mean_exchange, 
             time = times[time],
             amino_acid = split_sequence,
             unique_acid = paste0(1:length(split_sequence), split_sequence))
})



result_markov <- do.call(rbind, mean_exchange_markov)


ggplot(result_markov, aes(x = unique_acid, y = mean_exchange)) +
  geom_bar(stat = 'identity') + 
  scale_x_discrete(labels = split_sequence) +
  facet_wrap(~result$time) +
  ylim(0, 0.3)



# all_results_plot <- rbind(data.frame(result, sim = "from_package"), 
#                           data.frame(result_markov, sim = "markov"))
# 
# 
# ggplot(all_results_plot, aes(x = unique_acid, y = mean_exchange)) +
#   geom_bar(stat = 'identity') + 
#   scale_x_discrete(labels = split_sequence) +
#   facet_wrap(all_results_plot$sim~all_results_plot$time)







