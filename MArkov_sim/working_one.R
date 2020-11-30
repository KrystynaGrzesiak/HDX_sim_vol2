library(powerHDX)
library(ggplot2)
library(dplyr)
library(markovchain)
library(microbenchmark)

load("../theo_spectra_sim/sysdata.rda")


curves_using_package <- function(sequence, charge = NULL, protection_factor = 1,
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
    tryCatch({
      transition_probs = get_exchange_probabilities(kcHD, kcDH, deltaT, protection_factor)
      HD_matrices = get_HD_matrices(sequence, transition_probs,
                                    time_sequence, times,
                                    n_molecules)
      
      isotope_dists = lapply(1:length(times), function(ith_time) {
        observed_dist = get_observed_iso_dist(HD_matrices[[ith_time]], isotopic_probs, maxD)
        observed_peaks = matrix(0, maxD + maxND + 1, 2)
        DM = 1.00628
        observed_peaks[1, 1] = peptide_mass / charge + 1.007276
        observed_peaks[1, 2] = observed_dist[1]
        
        for (i in 2:(maxD + maxND + 1)) {
          observed_peaks[i, 1] = observed_peaks[i - 1, 1] + DM / charge
          observed_peaks[i, 2] = observed_dist[i]
        }
        data.frame(
          Exposure = times[ith_time],
          Mz = observed_peaks[, 1],
          Intensity = observed_peaks[, 2],
          PH = pH
        )
      })
      isotope_dists = do.call("rbind", isotope_dists)
    })
  }
  isotope_dists = rbind(data.frame(Exposure = 0,
                                   Mz = peptide_mass / charge + 1.007276,
                                   Intensity = isotopic_probs,
                                   PH = pH),
                        isotope_dists)
  isotope_dists[["Sequence"]] = paste0(sequence, collapse = "")
  if (length(unique(protection_factor)) == 1) {
    isotope_dists[["PF"]] = protection_factor[1]
  } else {
    isotope_dists[["PF"]] = paste(protection_factor,
                                  sep = ",", collapse = ",")
  }
  isotope_dists[["Charge"]] = charge
  isotope_dists[["Sim"]] = "pkg"
  isotope_dists = isotope_dists[isotope_dists[["Intensity"]] > min_probability, ]
  data.table::as.data.table(isotope_dists)
  get_deuteration_curve_single_spectrum(data.table::as.data.table(isotope_dists))
}



get_HD_matrices_using_markov = function(sequence, transition_probs, experiment_times,
                                        times_to_record, n_molecules = 100) {
  
  time_intervals <- cut(experiment_times, c(0, times_to_record), right = TRUE, include.lowest = TRUE)
  separated_times <- split(experiment_times, time_intervals)
  peptide_length = length(sequence)
  steps <- 0
  initial_state <- c(1, 0)

  hd_matrices <- lapply(1:length(times_to_record), function(i) {
    
    steps <- sum(lengths(separated_times)[1:i])
    HDmatrix <- sapply(1:peptide_length, function(amino_acid) {
      transition_matrix <- matrix(c(transition_probs[["HH"]][amino_acid], transition_probs[["HD"]][amino_acid],
                                    transition_probs[["DH"]][amino_acid], transition_probs[["DD"]][amino_acid]), nrow = 2, 
                                  byrow = TRUE)
      
      chain <- new('markovchain', states = c("0", "1"), 
                   transitionMatrix = transition_matrix, byrow = TRUE)  # 0 denotes hydrogen, 1 denotes deuterium
      sample(c(0, 1), n_molecules, replace = TRUE, prob = initial_state*chain^steps)
    })
    HDmatrix[, unique(c(1, 2, which(sequence == "P")))] = 0
    HDmatrix
  })
  hd_matrices
}


curves_using_markov_chain <- function(sequence, charge = NULL, protection_factor = 1,
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
    tryCatch({
      
      transition_probs = get_exchange_probabilities(kcHD, kcDH, deltaT, protection_factor)
      transition_probs[["HH"]] <- 1 - transition_probs[["HD"]]
      transition_probs[["DD"]] <- 1 - transition_probs[["DH"]]
      
      HD_matrices = get_HD_matrices_using_markov(sequence, transition_probs,
                                                 time_sequence, times,
                                                 n_molecules)
      
      isotope_dists = lapply(1:length(times), function(ith_time) {
        observed_dist = get_observed_iso_dist(HD_matrices[[ith_time]], isotopic_probs, maxD)
        observed_peaks = matrix(0, maxD + maxND + 1, 2)
        DM = 1.00628
        observed_peaks[1, 1] = peptide_mass / charge + 1.007276
        observed_peaks[1, 2] = observed_dist[1]
        
        for (i in 2:(maxD + maxND + 1)) {
          observed_peaks[i, 1] = observed_peaks[i - 1, 1] + DM / charge
          observed_peaks[i, 2] = observed_dist[i]
        }
        data.frame(
          Exposure = times[ith_time],
          Mz = observed_peaks[, 1],
          Intensity = observed_peaks[, 2],
          PH = pH
        )
      })
      isotope_dists = do.call("rbind", isotope_dists)
    })
  }
  isotope_dists = rbind(data.frame(Exposure = 0,
                                   Mz = peptide_mass / charge + 1.007276,
                                   Intensity = isotopic_probs,
                                   PH = pH),
                        isotope_dists)
  isotope_dists[["Sequence"]] = paste0(sequence, collapse = "")
  if (length(unique(protection_factor)) == 1) {
    isotope_dists[["PF"]] = protection_factor[1]
  } else {
    isotope_dists[["PF"]] = paste(protection_factor,
                                  sep = ",", collapse = ",")
  }
  isotope_dists[["Charge"]] = charge
  isotope_dists[["Sim"]] = "markov"
  isotope_dists = isotope_dists[isotope_dists[["Intensity"]] > min_probability, ]
  get_deuteration_curve_single_spectrum(data.table::as.data.table(isotope_dists))
}



