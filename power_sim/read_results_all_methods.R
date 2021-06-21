library(ggplot2)
library(dplyr)
library(stringr)

do_plot <- function(sim_results, pfs = unique(c(sim_results$State_1, sim_results$State_2))) {
  
  # New facet label names for dose variable
  pf1.labs <- paste("PF1 =", pfs)
  names(pf1.labs) <- as.character(pfs)
  
  # New facet label names for supp variable
  pf2.labs <- paste("PF2 =", pfs)
  names(pf2.labs) <- as.character(pfs)
  
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
          axis.title.x=element_blank(),
          legend.position="bottom") +
    theme(plot.title = element_text(hjust = 0.5))+
    geom_label(aes(label = formatC(Power, digits = 2, width = 1)), show.legend = FALSE)
  
  
  data_col <- sim_results %>% 
    filter(State_1 == State_2) %>% 
    filter(State_1 %in% pfs, State_2 %in% pfs)
  
  p + 
    geom_rect(data = data_col, col = "black", fill = "white", xmin = -Inf, xmax = Inf, 
              ymin = -Inf, ymax = Inf, alpha = 0)
}

read_summarize_files <- function(path) {
  files <- list.files(path = path, pattern = "\\.RDS$", full.names = TRUE) 
  do.call("rbind", lapply(files, readRDS)) %>% 
    filter(Transformation == "identity" | is.na(Transformation)) %>% 
    mutate(Test_id = paste0(Test, "_", Transformation)) %>%
    filter(Time != "categorical" | is.na(Time)) %>% 
    select(-Time)
}


sim_res <- read_summarize_files("./regularization/results")
sim_res_houde <- read_summarize_files("./houde_area/results")

tests <- c("Deuteros lm_identity", "Houde_NA","MEMHDX lmm_identity" ,
           "LASSO_knots_random_intercept_id_identity",
           "RIDGE_knots_random_intercept_id_exposure_identity",
           "RIDGE_knots_random_intercept_id_identity"
)

############ KOLORKI
col_plt1 <- c("thistle", "lightblue", "khaki", "hotpink", "dodgerblue1", "darkolivegreen2")

col_plt2 <- c("grey77", "grey47", "grey32", "dodgerblue1", "limegreen", "mediumpurple1")


col_plt3 <- c("red3", "firebrick2", "orangered", "dodgerblue1", "limegreen", "mediumpurple1")

col_plt3 <- c("firebrick4", "firebrick3", "firebrick1", "dodgerblue4", "dodgerblue3", "dodgerblue")

plot_rej_rate <- sim_res_houde %>% 
  filter(Sequence %in% intersect(Sequence, sim_res$Sequence)) %>% 
  rbind(sim_res) %>% 
  group_by(Test_id, State_1, State_2) %>% 
  summarise(Power = round(mean(Significant_difference, na.rm = TRUE), 2)) %>%
  ungroup() %>% 
  mutate(Test_id = factor(Test_id, levels = tests)) %>% 
  do_plot(pfs = c(20, 30, 40)) +
  # scale_fill_manual(name = "Test", 
  #                   values = col_plt3,
  #                   labels = c("HDX-Analyzer", "Houde","MEMHDX",  
  #                              "LASSO id", "RIDGE exposure id", "RIDGE id"))
  scale_fill_brewer(name = "Test", 
                    palette="Spectral", 
                    labels = c("HDX-Analyzer", "Houde","MEMHDX",  
                               "LASSO id", "RIDGE exposure id", "RIDGE id"))


#####################




plot_rej_rate <- sim_res_houde %>% 
  filter(Sequence %in% intersect(Sequence, sim_res$Sequence)) %>% 
  rbind(sim_res) %>% 
  group_by(Test_id, State_1, State_2) %>% 
  summarise(Power = round(mean(Significant_difference, na.rm = TRUE), 2)) %>%
  ungroup() %>% 
  mutate(Test_id = factor(Test_id, levels = tests)) %>% 
  do_plot(pfs = c(10, 15, 20, 30, 40))+
  scale_fill_brewer(name = "Test", 
                    palette="Spectral", 
                    labels = c("HDX-Analyzer", "Houde","MEMHDX",  
                               "LASSO id", "RIDGE exposure id", "RIDGE id"))





plot_rej_rate <- sim_res_houde %>% 
  filter(Sequence %in% intersect(Sequence, sim_res$Sequence)) %>% 
  rbind(sim_res) %>% 
  group_by(Test_id, State_1, State_2) %>% 
  summarise(Power = round(mean(Significant_difference, na.rm = TRUE), 2)) %>%
  ungroup() %>% 
  mutate(Test_id = factor(Test_id, levels = tests)) %>% 
  do_plot(pfs = c(90, 100, 110, 190, 200, 210)) +
  scale_fill_brewer(name = "Test", 
                    palette="Spectral", 
                    labels = c("HDX-Analyzer", "Houde","MEMHDX",  
                               "LASSO id", "RIDGE esposure id", "RIDGE id"))


plot_rej_rate <- sim_res_houde %>% 
  filter(Sequence %in% intersect(Sequence, sim_res$Sequence)) %>% 
  rbind(sim_res) %>% 
  group_by(Test_id, State_1, State_2) %>% 
  summarise(Power = round(mean(Significant_difference, na.rm = TRUE), 2)) %>%
  ungroup() %>% 
  mutate(Test_id = factor(Test_id, levels = tests)) %>% 
  do_plot(pfs = c(10, 50, 90, 300, 400)) +
  scale_fill_brewer(name = "Test", 
                    palette="Spectral", 
                    labels = c("HDX-Analyzer", "Houde","MEMHDX",  
                               "LASSO id", "RIDGE esposure id", "RIDGE id"))




#pfs all
c(10, 15, 20, 30, 40, 50, 90, 100, 110, 190, 200, 210, 300, 400)






