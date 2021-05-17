
library(ggplot2)
library(dplyr)
library(stringr)

files <- list.files(path = "./results", pattern = "\\.RDS$", full.names = TRUE) 
results <- do.call("rbind", lapply(files, readRDS))

results <- results %>% 
  filter(Transformation == "identity") %>% 
  mutate(Test_id = paste0(Test, "_", Transformation)) %>%
  filter(Time != "categorical" | is.na(Time)) %>% 
  select(-Time)

sim_results = results %>% 
  group_by(Test_id, State_1, State_2) %>% 
  summarise(Power = round(mean(Significant_difference, na.rm = TRUE), 2)) %>%
  ungroup()


select_results <- function(data, diff = 0.01) {
  
  tests_id <- data %>% 
    filter(State_1 == State_2) %>% 
    group_by(Test_id) %>% 
    summarise(m_error = mean(Power)) %>% 
    filter(m_error <= m_error[Test_id == "Deuteros lm_identity"] + diff) %>% 
    select(Test_id)
  
  data %>% 
    filter(Test_id %in% tests_id$Test_id)
}


p <- sim_results %>% 
  select_results(diff = 0.05) %>% 
  ggplot(aes(x = Test_id, y = Power, fill = Test_id)) +
  geom_col()+
  ylab("Rejection rate") +
  theme_minimal() +
  ylim(0, 1.15) +
  # geom_text(aes(label=Power), position=position_dodge(width=0.9), vjust=-0.25, size = 3) +
  facet_grid(State_1~State_2) +
  ggtitle("Rejection rate in pairwise testing") +
  theme(axis.text.x=element_blank(),
        axis.title.x=element_blank()) +
  theme(plot.title = element_text(hjust = 0.5))+
  geom_label(aes(label = formatC(Power, digits = 2, width = 4)))


data_col <- results %>% 
  mutate(Test_id = paste0(Test, "_", Transformation)) %>% 
  group_by(Test_id, State_1, State_2) %>% 
  summarise(Power = round(mean(Significant_difference, na.rm = TRUE), 2)) %>% 
  filter(State_1 == State_2)

p + 
  geom_rect(data = data_col, col = "black", fill = "white", xmin = -Inf, xmax = Inf, 
            ymin = -Inf, ymax = Inf, alpha = 0)




















pfs <- c(10, 15, 90, 100, 200)


p <- sim_results %>% 
  select_results(diff = 0.1) %>% 
  filter(State_1 %in% pfs, State_2 %in% pfs) %>% 
  ggplot(aes(x = Test_id, y = Power, fill = Test_id)) +
  geom_col()+
  ylab("Rejection rate") +
  theme_minimal() +
  ylim(0, 1.15) +
  # geom_text(aes(label=Power), position=position_dodge(width=0.9), vjust=-0.25, size = 3) +
  facet_grid(State_1~State_2) +
  ggtitle("Rejection rate in pairwise testing") +
  theme(axis.text.x=element_blank(),
        axis.title.x=element_blank()) +
  theme(plot.title = element_text(hjust = 0.5))+
  geom_label(aes(label = formatC(Power, digits = 2, width = 1)))


data_col <- results %>% 
  mutate(Test_id = paste0(Test, "_", Transformation)) %>% 
  filter(State_1 %in% pfs, State_2 %in% pfs) %>% 
  group_by(Test_id, State_1, State_2) %>% 
  summarise(Power = round(mean(Significant_difference, na.rm = TRUE), 2)) %>% 
  select_results(diff = 0.1) %>% 
  filter(State_1 == State_2)

p + 
  geom_rect(data = data_col, col = "black", fill = "white", xmin = -Inf, xmax = Inf, 
            ymin = -Inf, ymax = Inf, alpha = 0)










