
library(ggplot2)


files <- list.files(path = "G:/Matematyka/magisterka/Semiparametric regression/semiparametric_writeup/sim_res2", pattern = "\\.RDS$", full.names = TRUE) 
sim_results <- do.call(rbind, do.call("rbind", lapply(files, readRDS)))


results <- sim_results %>% 
  mutate(State_1 = str_remove(State_1, "_2"),
         State_2 = str_remove(State_2, "_2"))



results %>%
  filter(State_1 == State_2) %>% 
  group_by(Test, State_1, State_2) %>% 
  summarise(Error = round(mean(Significant_difference, na.rm = TRUE), 3)) %>%
  mutate(Pair = paste(State_1, "and", State_2)) %>%
  ggplot(aes(y = Error, x = Test, fill = Test)) +
  geom_bar(stat = 'identity', position = 'dodge') +
  geom_text(aes(label=Error), position=position_dodge(width=0.9), vjust=-0.25) +
  facet_grid( ~ Pair) +
  theme_minimal() +
  ylab("Estimated type I error") +
  theme(axis.text.x=element_blank(),
        axis.title.x=element_blank()) +
  ylim(0, 0.2)


results %>% 
  filter(State_1 != State_2) %>% 
  group_by(Test, State_1, State_2) %>% 
  summarise(Power = round(mean(Significant_difference, na.rm = TRUE), 2)) %>%
  mutate(Pair = paste(State_1, "and", State_2)) %>%
  ggplot(aes(y = Power, x = Test, fill = Test)) +
  geom_bar(stat = 'identity', position = 'dodge') +
  geom_text(aes(label=Power), position=position_dodge(width=0.9), vjust=-0.25) +
  facet_wrap( ~ Pair) +
  theme_minimal() +
  theme(axis.text.x=element_blank(),
        axis.title.x=element_blank()) +
  ylab("Estimated power") +
  ylim(0, 1.15)


results %>% 
  group_by(Test, State_1, State_2) %>% 
  filter(Test %in% c("Deuteros lm", "Semiparametric1", "Semiparametric4", "Semiparametric5")) %>% 
  summarise(Power = round(mean(Significant_difference, na.rm = TRUE), 2)) %>%
  ggplot(aes(x = Test, y = Power, fill = Test)) +
  geom_col()+
  ylab("Rejection rate") +
  theme_minimal() +
  ylim(0, 1.15) +
  geom_text(aes(label=Power), position=position_dodge(width=0.9), vjust=-0.25, size = 3) +
  facet_grid(State_1~State_2) +
  ggtitle("Rejection rate in pairwise testing") +
  theme(axis.text.x=element_blank(),
        axis.title.x=element_blank()) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette="Paired", labels = c("Deuteros lm", "Log interaction", "State",
                                                 "Interaction"))