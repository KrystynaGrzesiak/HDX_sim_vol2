library(doParallel)
library(powerHDX)

setwd("~/HDX_sim_vol2/theo_spectra_sim")
load("sysdata.rda")


look_for_big_vectors = function (sequence, charge = NULL, protection_factor = 1, times = c(60, 
                                                                                           600), pH = 7.5, temperature = 15, n_molecules = 100, time_step_const = 1, 
                                 if_corr = 0, min_probability = 1e-04) 
{
  sequence = strsplit(sequence, "")[[1]]
  if (length(protection_factor) == 1L) {
    protection_factor = rep(protection_factor, length(sequence))
  }
  if (is.null(charge)) {
    charge = sample(2:6, 1)
  }
  peptide_iso_dist = get_approx_isotopic_distribution(sequence, 
                                                      min_probability)
  peptide_mass = peptide_iso_dist[[1]]
  isotopic_probs = peptide_iso_dist[[2]]
  maxND = peptide_iso_dist[[3]]
  maxD = peptide_iso_dist[[4]]
  kcHD = get_exchange_rates(sequence, "HD", pH, temperature, 
                            "poly", if_corr)
  kcDH = get_exchange_rates(sequence, "DH", pH, temperature, 
                            "poly")
  kmax = max(max(kcDH), max(kcHD))
  deltaT = time_step_const/kmax
  time_sequence = seq(0, max(times), deltaT)
  if (time_sequence == 0 && length(time_sequence) == 1) {
    print("There is no deuteration before given time point.")
    isotope_dists = data.frame()
  }
  else {
    tryCatch({
      times_to_record = get_recording_times(time_sequence, 
                                            times)
      times_to_record = setdiff(times_to_record, 0)
      transition_probs = get_exchange_probabilities(kcHD, 
                                                    kcDH, deltaT, protection_factor)
      HD_matrices = get_HD_matrices(sequence, transition_probs, 
                                    time_sequence, times_to_record, n_molecules)
      isotope_dists = lapply(1:length(times_to_record), 
                             function(ith_time) {
                               observed_dist = get_observed_iso_dist(HD_matrices[[ith_time]], 
                                                                     isotopic_probs, maxD)
                               observed_peaks = matrix(0, maxD + maxND + 1, 
                                                       2)
                               DM = 1.00628
                               observed_peaks[1, 1] = peptide_mass/charge + 
                                 1.007276
                               observed_peaks[1, 2] = observed_dist[1]
                               for (i in 2:(maxD + maxND + 1)) {
                                 observed_peaks[i, 1] = observed_peaks[i - 
                                                                         1, 1] + DM/charge
                                 observed_peaks[i, 2] = observed_dist[i]
                               }
                               data.frame(Exposure = times[ith_time], Mz = observed_peaks[, 
                                                                                          1], Intensity = observed_peaks[, 2], PH = pH)
                             })
      isotope_dists = do.call("rbind", isotope_dists)
    }, print(paste(paste0(sequence), "pH:", pH, "PF : ", protection_factor, "charge: ", charge)))
  }
  isotope_dists = rbind(data.frame(Exposure = 0, Mz = peptide_mass/charge + 
                                     1.007276, Intensity = isotopic_probs, PH = pH), isotope_dists)
  isotope_dists[["Sequence"]] = paste0(sequence, collapse = "")
  if (length(unique(protection_factor)) == 1) {
    isotope_dists[["PF"]] = protection_factor[1]
  }
  else {
    isotope_dists[["PF"]] = paste(protection_factor, 
                                  sep = ",", collapse = ",")
  }
  isotope_dists[["Charge"]] = charge
  isotope_dists = isotope_dists[isotope_dists[["Intensity"]] > 
                                  min_probability, ]
  data.table::as.data.table(isotope_dists)
}



sim_theo_spectra = function(all_params, n_cpus, times) {
  all_params$cpu = sample(1:n_cpus, nrow(all_params), replace = TRUE)
  all_params = split(all_params, all_params$cpu)
  
  mclapply(all_params, function(df_per_cpu) {
    df_per_cpu$sequence = as.character(df_per_cpu$sequence)
    df_per_cpu$charge = as.numeric(as.character(df_per_cpu$charge))
    lapply(1:nrow(df_per_cpu), function(ith_row) {
      print(paste("Simulation", ith_row, "\n"))
      res = tryCatch(look_for_big_vectors(sequence = df_per_cpu[ith_row, "sequence"],
                                          charge = df_per_cpu[ith_row, "charge"],
                                          protection_factor = df_per_cpu[ith_row, "protection_factor"],
                                          times = times,
                                          pH = df_per_cpu[ith_row, "pH"],
                                          temperature = 15,
                                          n_molecules = 500,
                                          time_step_const = df_per_cpu[ith_row, "step"]),
                     error = function(e) {
                       print(e)
                       data.frame()})
    })
  }, mc.cores = n_cpus)
}

times = c(5, 10, 20, 30, 40, 50, 60, 100, 300, 500, 900, 1200, 1500, 1800,
          2100, 2400, 3600, 7200, 21600, 43200)


all_params = readRDS("./missing_params.RDS")
all_params$sequence = as.character(all_params$sequence)
all_params$charge = as.numeric(as.character(all_params$charge))

cores = detectCores()

sim_theo_spectra(all_params, cores, times)






