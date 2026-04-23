# Load libraries
library(cmdstanr)
library(bayesplot)
library(MCMCpack)
library(dplyr)
library(tidyverse)
library(posterior)
library(patchwork)
library(parallel)
library(latex2exp)
library(scoringutils)
library(qs2)

# run auxiliary functions
source("code/aux_fun_sim.R")

# load sample
sim1.1 <- qs_read("results/samples/normal/sim1_1_normal_fixed_var.qs2")
sim1.2 <- qs_read("results/samples/normal/sim1_2_normal_fixed_var.qs2")
sim1.3 <- qs_read("results/samples/normal/sim1_3_normal_fixed_var.qs2")
sim2.1 <- qs_read("results/samples/normal/sim2_1_normal_fixed_var.qs2")
sim2.2 <- qs_read("results/samples/normal/sim2_2_normal_fixed_var.qs2")
sim2.3 <- qs_read("results/samples/normal/sim2_3_normal_fixed_var.qs2")
sim3.1 <- qs_read("results/samples/normal/sim3_1_normal_fixed_var.qs2")
sim3.2 <- qs_read("results/samples/normal/sim3_2_normal_fixed_var.qs2")
sim3.3 <- qs_read("results/samples/normal/sim3_3_normal_fixed_var.qs2")
sim4.1 <- qs_read("results/samples/normal/sim4_1_normal_fixed_var.qs2")
sim4.2 <- qs_read("results/samples/normal/sim4_2_normal_fixed_var.qs2")
sim4.3 <- qs_read("results/samples/normal/sim4_3_normal_fixed_var.qs2")
sim5.1 <- qs_read("results/samples/normal/sim5_1_normal_fixed_var.qs2")

# combine simulations
sim_sces <- list(
  sim1.1, sim1.2, sim1.3,
  sim2.1, sim2.2, sim2.3,
  sim3.1, sim3.2, sim3.3,
  sim4.1, sim4.2, sim4.3,
  sim5.1
)

# Compute MSE
true_value <- 5
mses <- lapply(sim_sces, function(sim) {
  colMeans((sim$hattheta - true_value)^2)
})

biass <- lapply(sim_sces, function(sim) {
  colMeans(sim$hattheta) - true_value
})

# Compute BCI
bcis_95 <- lapply(sim_sces, function(sim) {
  compute_bci(sim$theta, 0.95)
})
bcis_90 <- lapply(sim_sces, function(sim) {
  compute_bci(sim$theta, 0.90)
})
bcis_80 <- lapply(sim_sces, function(sim) {
  compute_bci(sim$theta, 0.80)
})
bcis_50 <- lapply(sim_sces, function(sim) {
  compute_bci(sim$theta, 0.50)
})
# compute coverage
coverages_95 <- lapply(1:length(bcis_95), function(i) {
  compute_coverage(bcis_95[[i]], true_value)
})
coverages_90 <- lapply(1:length(bcis_90), function(i) {
  compute_coverage(bcis_90[[i]], true_value)
})
coverages_80 <- lapply(1:length(bcis_80), function(i) {
  compute_coverage(bcis_80[[i]], true_value)
})
coverages_50 <- lapply(1:length(bcis_50), function(i) {
  compute_coverage(bcis_50[[i]], true_value)
})

# compute WIS
quantile_level <- c(0.025, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.975)
wis_list <- lapply(sim_sces, function(sim) {
  compute_wis(sim$theta, quantile_level, true_value)
})


# plot BCI
plots_bci_95 <- lapply(1:length(bcis_95), function(i) {
  sce <- as.roman(ceiling(i/3))
  sce <- ifelse(i %% 3 == 1, paste0(sce,".I"), ifelse(i %% 3 == 2, paste0(sce,".II"), paste0(sce,".III")))
  title_left <- ggplot() +
    annotate("text", x = 0, y = 0,
             label = sce,
             angle = 90, size = 8) +
    theme_void()
  plot <- plot_bci(bcis_95[[i]], sim_sces[[i]], true_value, 0.95, coverages_95[[i]], sce, "Bernoulli")
  plot <- (title_left | plot) + plot_layout(widths = c(0.05, 1))
  if (i == length(bcis_95)) {
    plot <- plot & theme(legend.position = "bottom")
  }
  return(plot)
})
wrap_plots_95 <- Reduce(`/`, plots_bci_95)
ggsave("results/figures/normal/bci/bci_95_normal_fixed_var.pdf",
        wrap_plots_95, width = 21, height = 29, dpi = 300, limitsize = F)
plots_bci_90 <- lapply(1:length(bcis_90), function(i) {
  sce <- as.roman(ceiling(i/3))
  sce <- ifelse(i %% 3 == 1, paste0(sce,".I"), ifelse(i %% 3 == 2, paste0(sce,".II"), paste0(sce,".III")))
  title_left <- ggplot() +
    annotate("text", x = 0, y = 0,
             label = sce,
             angle = 90, size = 8) +
    theme_void()
  plot <- plot_bci(bcis_90[[i]], sim_sces[[i]], true_value, 0.90, coverages_90[[i]], sce, "Bernoulli")
  plot <- (title_left | plot) + plot_layout(widths = c(0.05, 1))
  if (i == length(bcis_95)) {
    plot <- plot & theme(legend.position = "bottom")
  }
  return(plot)
})
wrap_plots_90 <- Reduce(`/`, plots_bci_90)
ggsave("results/figures/normal/bci/bci_90_normal_fixed_var.pdf", 
        wrap_plots_90, width = 21, height = 29, dpi = 300, limitsize = F)
plots_bci_80 <- lapply(1:length(bcis_80), function(i) {
  sce <- as.roman(ceiling(i/3))
  sce <- ifelse(i %% 3 == 1, paste0(sce,".I"), ifelse(i %% 3 == 2, paste0(sce,".II"), paste0(sce,".III")))
  title_left <- ggplot() +
    annotate("text", x = 0, y = 0,
             label = sce,
             angle = 90, size = 8) +
    theme_void()
  plot <- plot_bci(bcis_80[[i]], sim_sces[[i]], true_value, 0.80, coverages_80[[i]], sce, "Bernoulli")
  plot <- (title_left | plot) + plot_layout(widths = c(0.05, 1))
  if (i == length(bcis_95)) {
    plot <- plot & theme(legend.position = "bottom")
  }
  return(plot)
})
wrap_plots_80 <- Reduce(`/`, plots_bci_80)
ggsave("results/figures/normal/bci/bci_80_normal_fixed_var.pdf",
        wrap_plots_80, width = 21, height = 29, dpi = 300, limitsize = F)
plots_bci_50 <- lapply(1:length(bcis_50), function(i) {
  sce <- as.roman(ceiling(i/3))
  sce <- ifelse(i %% 3 == 1, paste0(sce,".I"), ifelse(i %% 3 == 2, paste0(sce,".II"), paste0(sce,".III")))
  title_left <- ggplot() +
    annotate("text", x = 0, y = 0,
             label = sce,
             angle = 90, size = 8) +
    theme_void()
  plot <- plot_bci(bcis_50[[i]], sim_sces[[i]], true_value, 0.50, coverages_50[[i]], sce, "Bernoulli")
  plot <- (title_left | plot) + plot_layout(widths = c(0.05, 1))
  if (i == length(bcis_95)) {
    plot <- plot & theme(legend.position = "bottom")
  }
  return(plot)
})
wrap_plots_50 <- Reduce(`/`, plots_bci_50)
ggsave("results/figures/normal/bci/bci_50_normal_fixed_var.pdf",
        wrap_plots_50, width = 21, height = 29, dpi = 300, limitsize = F)


# ------------------------------------------------------------------------------
# Summary results
# ------------------------------------------------------------------------------

summary_table <- data.frame(scenario = unlist(lapply(1:length(sim_sces), 
                                                     function(i) {
                  sce <- as.roman(ceiling(i/3))
                  sce <- ifelse(i %% 3 == 1, paste0(sce,".I"), ifelse(i %% 3 == 2, paste0(sce,".II"), paste0(sce,".III")))
                  rep(sce,4)
                })
                ),
                model = rep(c("NPP", "NPP-SEQ", "ONPP", "ONPP-SEQ"),length(sim_sces)),
                bias = unlist(biass),
                mse = unlist(mses),
                average_wis = unlist(wis_list),
                average_len_95 = unlist(lapply(bcis_95, compute_avg_len)),
                average_len_90 = unlist(lapply(bcis_90, compute_avg_len)),
                average_len_80 = unlist(lapply(bcis_80, compute_avg_len)),
                average_len_50 = unlist(lapply(bcis_50, compute_avg_len)),
                coverage_95 = unlist(coverages_95),
                coverage_90 = unlist(coverages_90),
                coverage_80 = unlist(coverages_80),
                coverage_50 = unlist(coverages_50)
)
colnames(summary_table) <- c("Scenario", "Prior", "Bias", "MSE", "Avg_WIS", 
                              "Avg_Len_BCI95", "Avg_Len_BCI90", "Avg_Len_BCI80", "Avg_Len_BCI50",
                             "Cov_BCI95", "Cov_BCI90", "Cov_BCI80", "Cov_BCI50")
summary_table$Scenario <- ""
row.names(summary_table) <- NULL

print(xtable::xtable(summary_table %>% 
                       mutate(
                         MSE = formatC(MSE, format = "f", digits = 3),
                         Bias = formatC(Bias, format = "f", digits = 3),
                         Avg_WIS = formatC(Avg_WIS, format = "f", digits = 3)
                       )
), 
type = "latex", include.rownames = FALSE)
