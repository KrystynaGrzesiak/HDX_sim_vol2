
library(powerHDX)
library(dplyr)
library(ggplot2)
library(plyr)

files <- list.files(path = "G:/HDX_sim_vol2/theo_spectra_sim/results", pattern = "\\.RDS$", full.names = TRUE) 
sim_results <- do.call("rbind", lapply(files, readRDS))

params_from_file <- readRDS("all_params_new_pf.RDS") 