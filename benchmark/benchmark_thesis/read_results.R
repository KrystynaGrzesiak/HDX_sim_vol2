
library(ggplot2)
library(dplyr)
library(stringr)
library(tidyr)

files <- list.files(path = "./res_thesis", pattern = "\\.RDS$", full.names = TRUE) 
sim_results <- do.call("rbind", lapply(files, readRDS))


sim_results %>%
  mutate(markov_times = markov_times/1000000000,
         rcpp_times = rcpp_times/1000000000) %>% 
  rename(markov = markov_times,
         rcpp = rcpp_times) %>% 
  gather(algorithm, time, markov, rcpp) %>% 
  group_by(algorithm, intervals) %>% 
  mutate(m_time = mean(time)) %>% 
  ggplot(aes(x = factor(intervals), y = m_time, col = algorithm, group = algorithm)) +
  geom_point() +
  geom_line()
