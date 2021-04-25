
library(ggplot2)
library(dplyr)
library(stringr)

set.seed(10)
params <- readRDS("./all_params.RDS")
sequences <- sample(unique(params$sequence), 100, replace = FALSE)


files <- list.files(path = "./results", pattern = "\\.RDS$", full.names = TRUE) 
sim_results <- do.call("rbind", lapply(files, readRDS))


sim_results %>% 
  group_by(expr) %>% 
  summarise(mean = mean(time)/1000000000)
