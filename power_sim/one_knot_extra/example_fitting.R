
library(powerHDX)
library(xtable)
library(gridExtra)
library(grid)

truncated_lines <- function(x, kappa){
  (x - kappa)*(x > kappa)
}

set.seed(17)

times = c(5, 10, 20, 30, 40, 50, 60, 100, 300, 500, 900, 1200, 1500, 1800,
          2100, 2400, 3600, 7200, 21600, 43200)


spec1 <- simulate_theoretical_spectra("LVRKDLQN", protection_factor = 10, charge = 1:3, times = times)
spec2 <- simulate_theoretical_spectra("LVRKDLQN", protection_factor = 40, charge = 1:3, times = times)

spec <- rbind(spec1, spec2)

data <- get_noisy_deuteration_curves(spec, 
                                     compare_pairs = TRUE, 
                                     reference = "all", 
                                     n_replicates = 1)[[1]][[1]]


significance_level = 0.05
data[["id"]] <- paste0(data$Rep, data$Charge, data$Experimental_state)

States = unique(data$State)
K <- 3600
X <- truncated_lines(data$Exposure, K)


model = lmerTest::lmer(Mass ~ Exposure*State + (1|id) + X,
                       data = data,
                       REML = FALSE)

fitted <- fitted.values(model)

top_plot <- data %>% 
  mutate(fitted = fitted) %>% 
  rename("PF" = State) %>% 
  ggplot(aes(x = Exposure, y = Mass, col = PF)) +
  geom_point(alpha = 0.4) +
  geom_line(aes(x = Exposure, y = fitted, col = PF), size = 1) +
  theme_minimal() +
  scale_color_manual(values=c("chocolate2", "royalblue")) +
  theme(legend.position="bottom")


bottom_plot <- data %>% 
  mutate(fitted = fitted) %>% 
  rename("PF" = State) %>% 
  ggplot(aes(x = Exposure, y = Mass, col = PF)) +
  geom_point(alpha = 0.4) +
  geom_line(aes(x = Exposure, y = fitted, col = PF), size = 1) +
  theme_minimal() +
  xlim(0, 7200) +
  scale_color_manual(values=c("chocolate2", "royalblue")) +
  facet_wrap(~PF, labeller = label_both) +
  theme(legend.position='none')


plot_fitting <- grid.arrange(bottom_plot, top_plot, nrow = 2, 
                             top = textGrob("Sequence: LVRKDLQN. Knot at 3600 and random intercept among different curves", gp = gpar(fontsize = 15)))



### HDX Analyzer




significance_level = 0.05
data[["id"]] <- paste0(data$Rep, data$Charge, data$Experimental_state)

States = unique(data$State)

model = lm(Mass ~ Exposure*State,
                       data = data)

fitted <- fitted.values(model)

top_plot <- data %>% 
  mutate(fitted = fitted) %>% 
  rename("PF" = State) %>% 
  ggplot(aes(x = Exposure, y = Mass, col = PF)) +
  geom_point(alpha = 0.4) +
  geom_line(aes(x = Exposure, y = fitted, col = PF), size = 1) +
  theme_minimal() +
  scale_color_manual(values=c("chocolate2", "royalblue")) +
  theme(legend.position="bottom")


bottom_plot <- data %>% 
  mutate(fitted = fitted) %>% 
  rename("PF" = State) %>% 
  ggplot(aes(x = Exposure, y = Mass, col = PF)) +
  geom_point(alpha = 0.4) +
  geom_line(aes(x = Exposure, y = fitted, col = PF), size = 1) +
  theme_minimal() +
  xlim(0, 7200) +
  scale_color_manual(values=c("chocolate2", "royalblue")) +
  facet_wrap(~PF, labeller = label_both) +
  theme(legend.position='none')


plot_fitting <- grid.arrange(bottom_plot, top_plot, nrow = 2, 
                             top = textGrob("Sequence: LVRKDLQN. HDX - Analyzer", gp = gpar(fontsize = 15)))












