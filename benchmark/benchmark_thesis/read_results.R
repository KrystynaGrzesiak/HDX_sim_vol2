
library(ggplot2)
library(dplyr)
library(stringr)
library(tidyr)


files <- list.files(path = "./results", pattern = "\\.RDS$", full.names = TRUE) 
sim_results <- do.call("rbind", lapply(files, readRDS))
# 
# files <- list.files(path = "./res_more", pattern = "\\.RDS$", full.names = TRUE) 
# sim_results_m <- do.call("rbind", lapply(files, readRDS))


sim_results %>%
  mutate(markov_times = markov_times/1000000000,
         rcpp_times = rcpp_times/1000000000) %>% 
  rename(markov = markov_times,
         rcpp = rcpp_times) %>% 
  gather(algorithm, time, markov, rcpp) %>% 
  group_by(algorithm, intervals) %>% 
  summarise(m_time = mean(time)) %>% 
  filter(intervals != "(45,60]" | algorithm != "rcpp") %>% 
  ggplot(aes(x = factor(intervals), y = m_time, col = algorithm, group = algorithm)) +
  geom_point() +
  geom_line() +
  theme_minimal() +
  labs(x = "Length of the protein",
       y = "Mean time [seconds]") +
  ggtitle("Comparison of simulation times")+
  theme(plot.title = element_text(hjust = 0.5))





sim_results %>%
  mutate(markov_times = markov_times/1000000000,
         rcpp_times = rcpp_times/1000000000) %>% 
  rename(markov = markov_times,
         rcpp = rcpp_times) %>% 
  filter(markov < 0.2) %>% 
  gather(algorithm, time, markov, rcpp) %>% 
  group_by(algorithm, intervals) %>% 
  mutate(m_time = mean(time)) %>% 
  filter(intervals != "(45,60]" | algorithm != "rcpp") %>% 
  ggplot(aes(x = factor(intervals), y = time, fill = algorithm)) +
  geom_boxplot() +
  theme_minimal() +
  labs(x = "Length of the protein",
       y = "Time [seconds]") +
  ggtitle("Comparison of simulation times")+
  theme(plot.title = element_text(hjust = 0.5)) +
  facet_wrap(~algorithm, scales = "free_y") +
  theme(strip.background = element_blank(),
        strip.text.x = element_blank())


library(xtable)


sim_results %>% 
  select(-markov_times, -rcpp_times) %>% 
  unique() %>% 
  arrange(lengths) %>% 
  xtable()




sim_results %>%
  mutate(markov_times = markov_times/1000000000,
         rcpp_times = rcpp_times/1000000000) %>% 
  rename(markov = markov_times,
         rcpp = rcpp_times) %>% 
  filter(markov < 0.2) %>% 
  gather(algorithm, time, markov, rcpp) %>% 
  group_by(algorithm, intervals) %>% 
  summarise(m_time = mean(time)) %>% 
  filter(intervals != "(45,60]" | algorithm != "rcpp") %>% 
  xtable()

library(reshape2)

data_wide <- dcast(olddata_long, subject + sex ~ condition, value.var="measurement")
data_wide


sim_results %>%
  mutate(markov_times = markov_times/1000000000,
         rcpp_times = rcpp_times/1000000000) %>% 
  rename(markov = markov_times,
         rcpp = rcpp_times) %>% 
  filter(markov < 0.2) %>% 
  gather(algorithm, time, markov, rcpp) %>% 
  group_by(algorithm, intervals) %>% 
  summarise(sd = sd(time)/60) %>% 
  dcast(intervals ~ algorithm) %>% 
  xtable(digits = 5)










