
library(ggplot2)
library(dplyr)
library(stringr)

files <- list.files(path = "./results", pattern = "\\.RDS$", full.names = TRUE) 
sim_results <- do.call("rbind", lapply(files, readRDS))

sim_results = sim_results %>% 
  mutate(Test_id = paste0(Test, "_", Transformation))


p <- sim_results %>% 
  filter(Time != "categorical" | is.na(Time)) %>% 
  select(-Time) %>% 
  filter(!(Test %in% c("S4", "S1"))) %>% 
  group_by(Test_id, State_1, State_2) %>% 
  summarise(Power = round(mean(Significant_difference, na.rm = TRUE), 2)) %>%
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
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette="Paired")


data_col <- sim_results %>% 
  group_by(Test_id, State_1, State_2) %>% 
  summarise(Power = round(mean(Significant_difference, na.rm = TRUE), 2)) %>% 
  filter(State_1 == State_2)

p + 
  geom_rect(data = data_col, col = "black", fill = "white", xmin = -Inf, xmax = Inf, 
            ymin = -Inf, ymax = Inf, alpha = 0)



sim_results %>% 
  filter(Time != "categorical" | is.na(Time)) %>% 
  select(-Time) %>% 
  filter(!(Test %in% c("S4", "S1"))) %>% 
  filter(State_1 %in% c(10, 20, 30), State_2 %in% c(10, 20, 30)) %>% 
  group_by(Test_id, State_1, State_2) %>% 
  summarise(Power = round(mean(Significant_difference, na.rm = TRUE), 2)) %>%
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
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette="Paired") +
  geom_label(aes(label = formatC(Power, digits = 2, width = 4)))



