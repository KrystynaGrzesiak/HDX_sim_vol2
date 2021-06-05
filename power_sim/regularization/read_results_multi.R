
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


select_power <- function(data, diff = 0.01) {
  tests_id <- data %>% 
    filter(State_1 != State_2) %>% 
    group_by(Test_id) %>% 
    summarise(m_power = mean(Power)) %>% 
    filter(m_power >= m_power[Test_id == "Deuteros lm_identity"] + diff) %>% 
    select(Test_id)
  data %>% 
    filter(Test_id %in% c(tests_id$Test_id, "Deuteros lm_identity", "MEMHDX lmm_identity"))
}


p <- sim_results %>% 
  # select_results(diff = 0.01) %>% 
  # select_power(diff = 0.1) %>% 
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




















pfs <- unique(sim_results$State_1)


# New facet label names for dose variable
pf1.labs <- paste("PF1 =", pfs)
names(pf1.labs) <- as.character(pfs)

# New facet label names for supp variable
pf2.labs <- paste("PF2 =", pfs)
names(pf2.labs) <- as.character(pfs)


p <- sim_results %>% 
  select_results(diff = 0.01) %>% 
  select_power(diff = 0.1) %>% 
  filter(State_1 %in% pfs, State_2 %in% pfs) %>%
  ggplot(aes(x = Test_id, y = Power, fill = Test_id)) +
  geom_col()+
  ylab("Rejection rate") +
  theme_minimal() +
  ylim(0, 1.15) +
  # geom_text(aes(label=Power), position=position_dodge(width=0.9), vjust=-0.25, size = 3) +
  facet_grid(State_1~State_2, labeller = labeller(State_2 = pf1.labs, State_1 = pf2.labs)) +
  ggtitle("Rejection rate in pairwise testing") +
  theme(axis.text.x=element_blank(),
        axis.title.x=element_blank()) +
  theme(plot.title = element_text(hjust = 0.5))+
  geom_label(aes(label = formatC(Power, digits = 2, width = 1)))


data_col <- results %>% 
  mutate(Test_id = paste0(Test, "_", Transformation)) %>% 
  filter(Test_id %in% c("2100_id_identity", "2400_id_identity", "3600_id_identity", "Deuteros lm_identity")) %>% 
  filter(State_1 %in% pfs, State_2 %in% pfs) %>% 
  group_by(Test_id, State_1, State_2) %>% 
  summarise(Power = round(mean(Significant_difference, na.rm = TRUE), 2)) %>% 
  select_results(diff = 0.1) %>% 
  filter(State_1 == State_2)

p + 
  geom_rect(data = data_col, col = "black", fill = "white", xmin = -Inf, xmax = Inf, 
            ymin = -Inf, ymax = Inf, alpha = 0)














files_1 <- list.files(path = "./results", pattern = "\\.RDS$", full.names = TRUE)[1:248] 
results <- do.call("rbind", lapply(files, readRDS))

set.seed(10)

files_df <- data.frame(file = files_1,
                       group = sample(1:4, size = 248, replace = TRUE))

splitted_files <- split(files_df, f = files_df$group)

results <- do.call("rbind", lapply(splitted_files, function(files_by_group) {
  cbind(do.call("rbind", lapply(files_by_group$file, readRDS)),
        group = unique(files_by_group$group))
}))

results <- results %>% 
  filter(Transformation == "identity") %>% 
  mutate(Test_id = paste0(Test, "_", Transformation)) %>%
  filter(Time != "categorical" | is.na(Time)) %>% 
  select(-Time)

sim_results = results %>% 
  group_by(Test_id, State_1, State_2, group) %>% 
  summarise(Power = round(mean(Significant_difference, na.rm = TRUE), 2)) %>%
  ungroup()



do_plot <- function(sim_results, pf2.labs, pf1.labs) {
  p <- sim_results %>% 
    filter(State_1 %in% pfs, State_2 %in% pfs) %>%
    ggplot(aes(x = Test_id, y = Power, fill = Test_id)) +
    geom_col()+
    ylab("Rejection rate") +
    theme_minimal() +
    ylim(0, 1.15) +
    # geom_text(aes(label=Power), position=position_dodge(width=0.9), vjust=-0.25, size = 3) +
    facet_grid(State_1~State_2, labeller = labeller(State_2 = pf1.labs, State_1 = pf2.labs)) +
    ggtitle("Rejection rate in pairwise testing") +
    theme(axis.text.x=element_blank(),
          axis.title.x=element_blank()) +
    theme(plot.title = element_text(hjust = 0.5))+
    geom_label(aes(label = formatC(Power, digits = 2, width = 1)))
  
  
  data_col <- results %>% 
    mutate(Test_id = paste0(Test, "_", Transformation)) %>% 
    filter(Test_id %in% c("2100_id_identity", "2400_id_identity", "3600_id_identity", "Deuteros lm_identity")) %>% 
    filter(State_1 %in% pfs, State_2 %in% pfs) %>% 
    group_by(Test_id, State_1, State_2) %>% 
    summarise(Power = round(mean(Significant_difference, na.rm = TRUE), 2)) %>% 
    select_results(diff = 0.1) %>% 
    filter(State_1 == State_2)
  
  p + 
    geom_rect(data = data_col, col = "black", fill = "white", xmin = -Inf, xmax = Inf, 
              ymin = -Inf, ymax = Inf, alpha = 0)
}


pfs <- unique(sim_results$State_1)

# New facet label names for dose variable
pf1.labs <- paste("PF1 =", pfs)
names(pf1.labs) <- as.character(pfs)

# New facet label names for supp variable
pf2.labs <- paste("PF2 =", pfs)
names(pf2.labs) <- as.character(pfs)


p1 <- do_plot(sim_results[sim_results$group == 1, ], pf2.labs, pf1.labs)
p2 <- do_plot(sim_results[sim_results$group == 2, ], pf2.labs, pf1.labs)
p3 <- do_plot(sim_results[sim_results$group == 3, ], pf2.labs, pf1.labs)
p4 <- do_plot(sim_results[sim_results$group == 4, ], pf2.labs, pf1.labs)


library(gridExtra)
grid.arrange(p1, p2, p3, p4, ncol=2)












