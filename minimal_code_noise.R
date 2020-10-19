if (!require(powerHDX)) {
    devtools::install_github("mstaniak/powerHDX")
}

library(data.table)
library(powerHDX)
library(parallel)
data.table::setDTthreads(1)

time_constant = readRDS("./time_constant.RDS")
by_ph_seq = readRDS("./by_ph_seq.RDS")

possible_times = list(
    c(0, 10, 60, 300, 500, 1800, 7200, 21600, 43200),
    c(0, 5, 10, 60, 100, 300, 500, 1200, 1800, 2400, 3600, 7200, 21600, 432000)
)
time_varying_mass_deviations = list(
    c(10, 50, 50, 30, 30, 30, 20, 10, 10),
    c(10, 50, 50, 50, 50, 30, 20, 20, 20, 20, 20, 10, 10, 10)
)

possible_times = list(
    c(0, 10, 60, 300, 500),
    c(0, 5, 10, 60, 100, 300, 500, 1200)
)
time_varying_mass_deviations = list(
    c(10, 50, 50, 30, 30),
    c(10, 50, 50, 50, 50, 30, 20, 20)
)


get_power = function(spectra_list, n_cores) {
    mclapply(
        spectra_list, function(spectrum) {
            for (i in 1:nrow(time_constant)) {
                if (is.na(time_constant[i, "per_run_deviations"])) {
                    per_run_deviation = NULL
                } else {
                    per_run_deviation = time_constant[i, "per_run_deviations"]
                }
                spectrum = spectrum[Exposure %in% possible_times[[time_constant[i, "n_timepoints"]]], ]
                tests_results = tryCatch({
                    noisy_curves = get_noisy_deuteration_curves(spectrum, reference = 1,
                                                                n_runs = time_constant[i, "num_reps"],
                                                                mass_deviations = time_constant[i, "mass_deviations"],
                                                                per_run_deviations = per_run_deviation)
                    calculate_hdx_power(noisy_curves,
                                        list(deuteros, lme_model, auc_test, memhdx_model),
                                        0.05, summarized = FALSE)
                }, error = function(e) data.table::data.table())
                seq = as.character(unique(noisy_curves[[1]][[1]]$Sequence))
                n_runs = time_constant[i, "num_reps"]
                n_times = length(possible_times[[time_constant[i, "n_timepoints"]]])
                mass_deviations = time_constant[i, "mass_deviations"]
                saveRDS(tests_results, file = paste("./results/", seq, n_runs, n_times, mass_deviations, per_run_deviation, ".RDS", sep = "_"))
            }
        }, mc.cores = n_cores
    )

}

spectrum = data.table(spectrum)

spectrum = by_ph_seq[[1]]
i = 1
