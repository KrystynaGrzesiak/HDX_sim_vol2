
library(ggplot2)
library(dplyr)
library(stringr)

files <- list.files(path = "./results", pattern = "\\.RDS$", full.names = TRUE) 
sim_results <- do.call(rbind, do.call("rbind", lapply(files, readRDS)))


results <- sim_results %>% 
  mutate(State_1 = str_remove(State_1, "_2"),
         State_2 = str_remove(State_2, "_2"))


results %>% 
  group_by(Test, State_1, State_2) %>% 
  summarise(Power = round(mean(Significant_difference, na.rm = TRUE), 2)) %>%
  ggplot(aes(x = Test, y = Power, fill = Test)) +
  geom_col()+
  ylab("Rejection rate") +
  theme_minimal() +
  ylim(0, 1.15) +
  geom_text(aes(label=Power), position=position_dodge(width=0.9), vjust=-0.25, size = 3) +
  facet_grid(as.numeric(State_1)~as.numeric(State_2)) +
  ggtitle("Rejection rate in pairwise testing") +
  theme(axis.text.x=element_blank(),
        axis.title.x=element_blank()) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette="Paired", labels = c("Deuteros lm", "Interaction"))
