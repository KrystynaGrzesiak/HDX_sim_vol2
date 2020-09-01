library(powerHDX)

load("sysdata.rda")

times = c(5, 10, 20, 30, 40, 50, 60, 100, 300, 500, 900, 1200, 1500, 1800,
          2100, 2400, 3600, 7200, 21600, 43200)


params = readRDS("./all_params.RDS")[1:10, ]
params$sequence = as.character(params$sequence)
params$charge = as.numeric(as.character(params$charge))


spectra_list = lapply(1:nrow(params), function(ith_row) {
    simulate_theoretical_spectra(sequence = params[ith_row, "sequence"],
                                 charge = params[ith_row, "charge"],
                                 protection_factor = params[ith_row, "protection_factor"],
                                 times = times,
                                 pH = params[ith_row, "pH"],
                                 temperature = 15,
                                 n_molecules = 100,
                                 time_step_const = params[ith_row, "step"])
})

s = rbindlist(spectra_list)
spectra_list = split(s, s$Sequence)
spectrum = spectra_list[[1]]
spectrum = spectrum[Exposure %in% possible_times[[time_constant[i, "n_timepoints"]]], ]

noisy_curves = get_noisy_deuteration_curves(spectrum, reference = 1, compare_pairs = FALSE)



power = lapply(spectra_list, function(spectrum) {
    for (i in 1:nrow(time_constant)) {
        if (is.na(time_constant[i, "per_run_deviations"])) {
            per_run_deviation = NULL
        } else {
            per_run_deviation = time_constant[i, "per_run_deviations"]
        }
        spectrum = spectrum[Exposure %in% possible_times[[time_constant[i, "n_timepoints"]]], ]
        tests_results = tryCatch({
            noisy_curves = get_noisy_deuteration_curves(spectrum, reference = 1, compare_pairs = FALSE,
                                                        n_runs = time_constant[i, "num_reps"],
                                                        mass_deviations = time_constant[i, "mass_deviations"],
                                                        per_run_deviations = per_run_deviation)
            calculate_hdx_power(noisy_curves,
                                list(deuteros),
                                0.05, summarized = FALSE)
        })
    }
})
