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

# run auxiliary functions
source("code/aux_fun_sim.R")

# load sample
load("results/samples/bin/sim1_1_bin.RData")
load("results/samples/bin/sim1_2_bin.RData")
load("results/samples/bin/sim1_3_bin.RData")
load("results/samples/bin/sim2_1_bin.RData")
load("results/samples/bin/sim2_2_bin.RData")
load("results/samples/bin/sim2_3_bin.RData")
load("results/samples/bin/sim3_1_bin.RData")
load("results/samples/bin/sim3_2_bin.RData")
load("results/samples/bin/sim3_3_bin.RData")
load("results/samples/bin/sim4_1_bin.RData")
load("results/samples/bin/sim4_2_bin.RData")
load("results/samples/bin/sim4_3_bin.RData")

# combine simulations
sim_sces <- list(
  sim1.1, sim1.2, sim1.3,
  sim2.1, sim2.2, sim2.3,
  sim3.1, sim3.2, sim3.3,
  sim4.1, sim4.2, sim4.3
)

# Compute MSE
true_value <- 0.5
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
quantile_level <- c(0.05, 0.1, 0.2, 0.5, 0.8, 0.9, 0.95)
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
ggsave("results/figures/bin/bci/bci_95_bin.pdf", wrap_plots_95, width = 21, height = 29, dpi = 300, limitsize = F)
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
ggsave("results/figures/bin/bci/bci_90_bin.pdf", wrap_plots_90, width = 21, height = 29, dpi = 300, limitsize = F)
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
ggsave("results/figures/bin/bci/bci_80_bin.pdf", wrap_plots_80, width = 21, height = 29, dpi = 300, limitsize = F)
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
ggsave("results/figures/bin/bci/bci_50_bin.pdf", wrap_plots_50, width = 21, height = 29, dpi = 300, limitsize = F)


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

# ------------------------------------------------------------------------------
# Plot boxplots
# ------------------------------------------------------------------------------

plots_bin <- lapply(1:length(sim_sces), plot_sce_bin)
lapply(1:4, function(i) {
  sce_plots <- lapply(1:3, function(j) {
    plots_bin[[(i-1)*3 + j]]
  })
  combined_plot <- Reduce(`/`, sce_plots)
  sce <- as.roman(i)
  ggsave(paste0("results/figures/bin/boxplot_sce_", sce, "_bin.pdf"), combined_plot, width = 10, height = 12, dpi = 300)
})


# ------------------------------------------------------------------------------
# Count how many times NPP and NPP-SEQ get the order right
# ------------------------------------------------------------------------------

order_list <- list(c(1,3),
                   c(1,3),
                   c(1,3),
                   c(1,2),
                   c(1,2),
                   c(1,2),
                   c(3,1),
                   c(3,1),
                   c(3,1)
)

order_counts <- lapply(1:length(order_list), function(i) {
  sim <- sim_sces[[i]]
  min_max <- order_list[[i]]
  delta_npp <- lapply(sim$delta, function(x) x$delta_npp)
  delta_nppseq <- lapply(sim$delta, function(x) x$delta_nppseq)
  count_npp <- lapply(delta_npp, function(x) {
    (apply(x, 1, which.min) == min_max[1]) &
      (apply(x, 1, which.max) == min_max[2])
  })
  avg_count_npp <- unlist(count_npp) %>% mean()
  count_nppseq <- lapply(delta_nppseq, function(x) {
    (apply(x, 1, which.min) == min_max[1]) &
      (apply(x, 1, which.max) == min_max[2])
  })
  avg_count_nppseq <- unlist(count_nppseq) %>% mean()
  return(list(
    npp = avg_count_npp,
    nppseq = avg_count_nppseq
  ))
})
order_counts_df <- data.frame(
  scenario = unlist(lapply(1:length(order_list), 
                           function(i) {
                             sce <- as.roman(ceiling(i/3))
                             sce <- ifelse(i %% 3 == 1, paste0(sce,".I"), ifelse(i %% 3 == 2, paste0(sce,".II"), paste0(sce,".III")))
                             rep(sce,2)
                           })
  ),
  model = rep(c("NPP", "NPP-SEQ"), length(order_list)),
  order_count = unlist(lapply(order_counts, function(x) c(x$npp, x$nppseq)))
)    
print(xtable::xtable(order_counts_df), 
      type = "latex", include.rownames = FALSE)

# ------------------------------------------------------------------------------
# Compute how badly ONPP and ONPP-SEQ estimate the order
# ------------------------------------------------------------------------------

order_error_list <- lapply(1:12, function(i) {
  sim <- sim_sces[[i]]
  delta_onpp <- lapply(sim$delta, function(x) x$delta_onpp)
  delta_onppseq <- lapply(sim$delta, function(x) x$delta_onppseq)
  avg_dist1_onpp <- unlist(lapply(delta_onpp, function(x) {
    apply(x, 1, function(y) {
      y[2] - y[1]
    })
  })) %>% mean()
  avg_dist2_onpp <- unlist(lapply(delta_onpp, function(x) {
    apply(x, 1, function(y) {
      y[3] - y[2]
    })
  })) %>% mean()
  avg_dist1_onppseq <- unlist(lapply(delta_onppseq, function(x) {
    apply(x, 1, function(y) {
      y[2] - y[1]
    })
  })) %>% mean()
  avg_dist2_onppseq <- unlist(lapply(delta_onppseq, function(x) {
    apply(x, 1, function(y) {
      y[3] - y[2]
    })
  })) %>% mean()
  return(list(
    onpp = c(avg_dist1_onpp, avg_dist2_onpp),
    onppseq = c(avg_dist1_onppseq, avg_dist2_onppseq)
  ))
})

order_error_df <- data.frame(
  scenario = unlist(lapply(1:12, 
                           function(i) {
                             sce <- as.roman(ceiling(i/3))
                             sce <- ifelse(i %% 3 == 1, paste0(sce,".I"), ifelse(i %% 3 == 2, paste0(sce,".II"), paste0(sce,".III")))
                             rep(sce,2)
                           })
  ),
  model = rep(c("ONPP", "ONPP-SEQ"), 12),
  dist_med_min = unlist(lapply(order_error_list, function(x) c(x$onpp[1], x$onppseq[1]))),
  dist_max_med = unlist(lapply(order_error_list, function(x) c(x$onpp[2], x$onppseq[2])))
)  
print(xtable::xtable(order_error_df %>%
                       mutate(
                         dist_med_min = formatC(dist_med_min, format = "f", digits = 3),
                         dist_max_med = formatC(dist_max_med, format = "f", digits = 3)
                       )),
      type = "latex", include.rownames = FALSE)


# ------------------------------------------------------------------------------

plot_delta1_3.1 <- ggplot() +
  geom_density(aes(x = sim3.1$delta[[1]]$delta_npp[,1], color = "NPP"), linewidth = 0.8) +
  geom_density(aes(x = sim3.1$delta[[1]]$delta_nppseq[,1], color = "NPP_SEQ"), linewidth = 0.8) +
  geom_density(aes(x = sim3.1$delta[[1]]$delta_onpp[,1], color = "ONPP"), linewidth = 0.8) +
  geom_density(aes(x = sim3.1$delta[[1]]$delta_onppseq[,1], color = "ONPP_SEQ"), linewidth = 0.8) +
  scale_color_manual(name = NULL, 
                     values = c("NPP" = RColorBrewer::brewer.pal(4, "Set1")[1], 
                                "NPP_SEQ" = RColorBrewer::brewer.pal(4, "Set1")[2],
                                "ONPP" = RColorBrewer::brewer.pal(4, "Set1")[3],
                                "ONPP_SEQ" = RColorBrewer::brewer.pal(4, "Set1")[4])) +
  labs(x = expression(delta[1]),
       y = NULL) +
  ggtitle("Scenario III.I")
plot_delta1_3.1
plot_delta2_3.1 <- ggplot() +
  geom_density(aes(x = sim3.1$delta[[1]]$delta_npp[,2], color = "NPP"), linewidth = 0.8) +
  geom_density(aes(x = sim3.1$delta[[1]]$delta_nppseq[,2], color = "NPP_SEQ"), linewidth = 0.8) +
  geom_density(aes(x = sim3.1$delta[[1]]$delta_onpp[,2], color = "ONPP"), linewidth = 0.8) +
  geom_density(aes(x = sim3.1$delta[[1]]$delta_onppseq[,2], color = "ONPP_SEQ"), linewidth = 0.8) +
  scale_color_manual(name = NULL, 
                     values = c("NPP" = RColorBrewer::brewer.pal(4, "Set1")[1], 
                                "NPP_SEQ" = RColorBrewer::brewer.pal(4, "Set1")[2],
                                "ONPP" = RColorBrewer::brewer.pal(4, "Set1")[3],
                                "ONPP_SEQ" = RColorBrewer::brewer.pal(4, "Set1")[4])) +
  labs(x = expression(delta[2]),
       y = NULL)
plot_delta2_3.1
plot_delta3_3.1 <- ggplot() +
  geom_density(aes(x = sim3.1$delta[[1]]$delta_npp[,3], color = "NPP"), linewidth = 0.8) +
  geom_density(aes(x = sim3.1$delta[[1]]$delta_nppseq[,3], color = "NPP_SEQ"), linewidth = 0.8) +
  geom_density(aes(x = sim3.1$delta[[1]]$delta_onpp[,3], color = "ONPP"), linewidth = 0.8) +
  geom_density(aes(x = sim3.1$delta[[1]]$delta_onppseq[,3], color = "ONPP_SEQ"), linewidth = 0.8) +
  scale_color_manual(name = NULL, 
                     values = c("NPP" = RColorBrewer::brewer.pal(4, "Set1")[1], 
                                "NPP_SEQ" = RColorBrewer::brewer.pal(4, "Set1")[2],
                                "ONPP" = RColorBrewer::brewer.pal(4, "Set1")[3],
                                "ONPP_SEQ" = RColorBrewer::brewer.pal(4, "Set1")[4])) +
  labs(x = expression(delta[3]),
       y = NULL)
plot_delta3_3.1
combined_delta_plots_3.1 <- (plot_delta1_3.1 + plot_delta2_3.1 + plot_delta3_3.1) + 
  plot_layout(ncol = 3, guides = "collect")
combined_delta_plots_3.1

plot_delta1_3.3 <- ggplot() +
  geom_density(aes(x = sim3.3$delta[[1]]$delta_npp[,1], color = "NPP"), linewidth = 0.8) +
  geom_density(aes(x = sim3.3$delta[[1]]$delta_nppseq[,1], color = "NPP_SEQ"), linewidth = 0.8) +
  geom_density(aes(x = sim3.3$delta[[1]]$delta_onpp[,1], color = "ONPP"), linewidth = 0.8) +
  geom_density(aes(x = sim3.3$delta[[1]]$delta_onppseq[,1], color = "ONPP_SEQ"), linewidth = 0.8) +
  scale_color_manual(name = NULL, 
                     values = c("NPP" = RColorBrewer::brewer.pal(4, "Set1")[1], 
                                "NPP_SEQ" = RColorBrewer::brewer.pal(4, "Set1")[2],
                                "ONPP" = RColorBrewer::brewer.pal(4, "Set1")[3],
                                "ONPP_SEQ" = RColorBrewer::brewer.pal(4, "Set1")[4])) +
  labs(x = expression(delta[1]),
       y = NULL) +
  ggtitle("Scenario III.III")
plot_delta1_3.3
plot_delta2_3.3 <- ggplot() +
  geom_density(aes(x = sim3.3$delta[[1]]$delta_npp[,2], color = "NPP"), linewidth = 0.8) +
  geom_density(aes(x = sim3.3$delta[[1]]$delta_nppseq[,2], color = "NPP_SEQ"), linewidth = 0.8) +
  geom_density(aes(x = sim3.3$delta[[1]]$delta_onpp[,2], color = "ONPP"), linewidth = 0.8) +
  geom_density(aes(x = sim3.3$delta[[1]]$delta_onppseq[,2], color = "ONPP_SEQ"), linewidth = 0.8) +
  scale_color_manual(name = NULL, 
                     values = c("NPP" = RColorBrewer::brewer.pal(4, "Set1")[1], 
                                "NPP_SEQ" = RColorBrewer::brewer.pal(4, "Set1")[2],
                                "ONPP" = RColorBrewer::brewer.pal(4, "Set1")[3],
                                "ONPP_SEQ" = RColorBrewer::brewer.pal(4, "Set1")[4])) +
  labs(x = expression(delta[2]),
       y = NULL)
plot_delta2_3.3
plot_delta3_3.3 <- ggplot() +
  geom_density(aes(x = sim3.3$delta[[1]]$delta_npp[,3], color = "NPP"), linewidth = 0.8) +
  geom_density(aes(x = sim3.3$delta[[1]]$delta_nppseq[,3], color = "NPP_SEQ"), linewidth = 0.8) +
  geom_density(aes(x = sim3.3$delta[[1]]$delta_onpp[,3], color = "ONPP"), linewidth = 0.8) +
  geom_density(aes(x = sim3.3$delta[[1]]$delta_onppseq[,3], color = "ONPP_SEQ"), linewidth = 0.8) +
  scale_color_manual(name = NULL, 
                     values = c("NPP" = RColorBrewer::brewer.pal(4, "Set1")[1], 
                                "NPP_SEQ" = RColorBrewer::brewer.pal(4, "Set1")[2],
                                "ONPP" = RColorBrewer::brewer.pal(4, "Set1")[3],
                                "ONPP_SEQ" = RColorBrewer::brewer.pal(4, "Set1")[4])) +
  labs(x = expression(delta[3]),
       y = NULL)
plot_delta3_3.3
combined_delta_plots_3.3 <- (plot_delta1_3.3 + plot_delta2_3.3 + plot_delta3_3.3) + 
  plot_layout(ncol = 3, guides = "collect")
combined_delta_plots_3.3



plot_npp_3.1 <- ggplot() +
  geom_density(aes(x = sim3.1$delta[[1]]$delta_npp[,1], color = "delta1"), linewidth = 0.8) +
  geom_density(aes(x = sim3.1$delta[[1]]$delta_npp[,2], color = "delta2"), linewidth = 0.8) +
  geom_density(aes(x = sim3.1$delta[[1]]$delta_npp[,3], color = "delta3"), linewidth = 0.8) +
  scale_color_manual(name = NULL, 
                     values = c("delta1" = RColorBrewer::brewer.pal(3, "Set1")[1], 
                                "delta2" = RColorBrewer::brewer.pal(3, "Set1")[2],
                                "delta3" = RColorBrewer::brewer.pal(3, "Set1")[3]),
                     labels = c(expression(delta[1]), 
                                expression(delta[2]),
                                expression(delta[3]))) +
  labs(x = expression(delta),
       y = NULL) +
  ggtitle(label = "Scenario III.I",
          subtitle = "NPP")

plot_nppseq_3.1 <- ggplot() +
  geom_density(aes(x = sim3.1$delta[[1]]$delta_nppseq[,1], color = "delta1"), linewidth = 0.8) +
  geom_density(aes(x = sim3.1$delta[[1]]$delta_nppseq[,2], color = "delta2"), linewidth = 0.8) +
  geom_density(aes(x = sim3.1$delta[[1]]$delta_nppseq[,3], color = "delta3"), linewidth = 0.8) +
  scale_color_manual(name = NULL, 
                     values = c("delta1" = RColorBrewer::brewer.pal(3, "Set1")[1], 
                                "delta2" = RColorBrewer::brewer.pal(3, "Set1")[2],
                                "delta3" = RColorBrewer::brewer.pal(3, "Set1")[3]),
                     labels = c(expression(delta[1]), 
                                expression(delta[2]),
                                expression(delta[3]))) +
  labs(x = expression(delta),
       y = NULL) +
  ggtitle(label = NULL,
          subtitle =  "NPP-SEQ")

plot_onpp_3.1 <- ggplot() +
  geom_density(aes(x = sim3.1$delta[[1]]$delta_onpp[,1], color = "delta1"), linewidth = 0.8) +
  geom_density(aes(x = sim3.1$delta[[1]]$delta_onpp[,2], color = "delta2"), linewidth = 0.8) +
  geom_density(aes(x = sim3.1$delta[[1]]$delta_onpp[,3], color = "delta3"), linewidth = 0.8) +
  scale_color_manual(name = NULL, 
                     values = c("delta1" = RColorBrewer::brewer.pal(3, "Set1")[1], 
                                "delta2" = RColorBrewer::brewer.pal(3, "Set1")[2],
                                "delta3" = RColorBrewer::brewer.pal(3, "Set1")[3]),
                     labels = c(expression(delta[1]), 
                                expression(delta[2]),
                                expression(delta[3]))) +
  labs(x = expression(delta),
       y = NULL) +
  ggtitle(label = NULL,
          subtitle = "ONPP")

plot_onppseq_3.1 <- ggplot() +
  geom_density(aes(x = sim3.1$delta[[1]]$delta_onppseq[,1], color = "delta1"), linewidth = 0.8) +
  geom_density(aes(x = sim3.1$delta[[1]]$delta_onppseq[,2], color = "delta2"), linewidth = 0.8) +
  geom_density(aes(x = sim3.1$delta[[1]]$delta_onppseq[,3], color = "delta3"), linewidth = 0.8) +
  scale_color_manual(name = NULL, 
                     values = c("delta1" = RColorBrewer::brewer.pal(3, "Set1")[1], 
                                "delta2" = RColorBrewer::brewer.pal(3, "Set1")[2],
                                "delta3" = RColorBrewer::brewer.pal(3, "Set1")[3]),
                     labels = c(expression(delta[1]), 
                                expression(delta[2]),
                                expression(delta[3]))) +
  labs(x = expression(delta),
       y = NULL) +
  ggtitle(label = NULL,
          subtitle = "ONPP-SEQ")
combined_models_3.1 <- (plot_npp_3.1 + plot_nppseq_3.1 + plot_onpp_3.1 + plot_onppseq_3.1) + 
  plot_layout(ncol = 4, guides = "collect")
combined_models_3.1

plot_npp_3.3 <- ggplot() +
  geom_density(aes(x = sim3.3$delta[[1]]$delta_npp[,1], color = "delta1"), linewidth = 0.8) +
  geom_density(aes(x = sim3.3$delta[[1]]$delta_npp[,2], color = "delta2"), linewidth = 0.8) +
  geom_density(aes(x = sim3.3$delta[[1]]$delta_npp[,3], color = "delta3"), linewidth = 0.8) +
  scale_color_manual(name = NULL, 
                     values = c("delta1" = RColorBrewer::brewer.pal(3, "Set1")[1], 
                                "delta2" = RColorBrewer::brewer.pal(3, "Set1")[2],
                                "delta3" = RColorBrewer::brewer.pal(3, "Set1")[3]),
                     labels = c(expression(delta[1]), 
                                expression(delta[2]),
                                expression(delta[3]))) +
  labs(x = expression(delta),
       y = NULL) +
  ggtitle(label = "Scenario III.III",
          subtitle = "NPP")
plot_nppseq_3.3 <- ggplot() +
  geom_density(aes(x = sim3.3$delta[[1]]$delta_nppseq[,1], color = "delta1"), linewidth = 0.8) +
  geom_density(aes(x = sim3.3$delta[[1]]$delta_nppseq[,2], color = "delta2"), linewidth = 0.8) +
  geom_density(aes(x = sim3.3$delta[[1]]$delta_nppseq[,3], color = "delta3"), linewidth = 0.8) +
  scale_color_manual(name = NULL, 
                     values = c("delta1" = RColorBrewer::brewer.pal(3, "Set1")[1], 
                                "delta2" = RColorBrewer::brewer.pal(3, "Set1")[2],
                                "delta3" = RColorBrewer::brewer.pal(3, "Set1")[3]),
                     labels = c(expression(delta[1]), 
                                expression(delta[2]),
                                expression(delta[3]))) +
  labs(x = expression(delta),
       y = NULL) +
  ggtitle(label = NULL,
          subtitle =  "NPP-SEQ")
plot_onpp_3.3 <- ggplot() +
  geom_density(aes(x = sim3.3$delta[[1]]$delta_onpp[,1], color = "delta1"), linewidth = 0.8) +
  geom_density(aes(x = sim3.3$delta[[1]]$delta_onpp[,2], color = "delta2"), linewidth = 0.8) +
  geom_density(aes(x = sim3.3$delta[[1]]$delta_onpp[,3], color = "delta3"), linewidth = 0.8) +
  scale_color_manual(name = NULL, 
                     values = c("delta1" = RColorBrewer::brewer.pal(3, "Set1")[1], 
                                "delta2" = RColorBrewer::brewer.pal(3, "Set1")[2],
                                "delta3" = RColorBrewer::brewer.pal(3, "Set1")[3]),
                     labels = c(expression(delta[1]), 
                                expression(delta[2]),
                                expression(delta[3]))) +
  labs(x = expression(delta),
       y = NULL) +
  ggtitle(label = NULL,
          subtitle = "ONPP")
plot_onppseq_3.3 <- ggplot() +
  geom_density(aes(x = sim3.3$delta[[1]]$delta_onppseq[,1], color = "delta1"), linewidth = 0.8) +
  geom_density(aes(x = sim3.3$delta[[1]]$delta_onppseq[,2], color = "delta2"), linewidth = 0.8) +
  geom_density(aes(x = sim3.3$delta[[1]]$delta_onppseq[,3], color = "delta3"), linewidth = 0.8) +
  scale_color_manual(name = NULL, 
                     values = c("delta1" = RColorBrewer::brewer.pal(3, "Set1")[1], 
                                "delta2" = RColorBrewer::brewer.pal(3, "Set1")[2],
                                "delta3" = RColorBrewer::brewer.pal(3, "Set1")[3]),
                     labels = c(expression(delta[1]), 
                                expression(delta[2]),
                                expression(delta[3]))) +
  labs(x = expression(delta),
       y = NULL) +
  ggtitle(label = NULL,
          subtitle = "ONPP-SEQ")
combined_models_3.3 <- (plot_npp_3.3 + plot_nppseq_3.3 + plot_onpp_3.3 + plot_onppseq_3.3) + 
  plot_layout(ncol = 4, guides = "collect")
combined_models_3.3

plot_hatdelta_npp_3.1 <- ggplot() +
  geom_density(aes(x = sim3.1$hatdelta[[1]]$hatdelta_npp, color = "delta1"), linewidth = 0.8) +
  geom_density(aes(x = sim3.1$hatdelta[[2]]$hatdelta_npp, color = "delta2"), linewidth = 0.8) +
  geom_density(aes(x = sim3.1$hatdelta[[3]]$hatdelta_npp, color = "delta3"), linewidth = 0.8) +
  scale_color_manual(name = NULL, 
                     values = c("delta1" = RColorBrewer::brewer.pal(3, "Set1")[1], 
                                "delta2" = RColorBrewer::brewer.pal(3, "Set1")[2],
                                "delta3" = RColorBrewer::brewer.pal(3, "Set1")[3]),
                     labels = c(expression(delta[1]), 
                                expression(delta[2]),
                                expression(delta[3]))) +
  labs(x = expression(bar(delta)),
       y = NULL) +
  ggtitle(label = "Scenario III.I",
          subtitle = "NPP")
plot_hatdelta_nppseq_3.1 <- ggplot() +
  geom_density(aes(x = sim3.1$hatdelta[[1]]$hatdelta_nppseq, color = "delta1"), linewidth = 0.8) +
  geom_density(aes(x = sim3.1$hatdelta[[2]]$hatdelta_nppseq, color = "delta2"), linewidth = 0.8) +
  geom_density(aes(x = sim3.1$hatdelta[[3]]$hatdelta_nppseq, color = "delta3"), linewidth = 0.8) +
  scale_color_manual(name = NULL, 
                     values = c("delta1" = RColorBrewer::brewer.pal(3, "Set1")[1], 
                                "delta2" = RColorBrewer::brewer.pal(3, "Set1")[2],
                                "delta3" = RColorBrewer::brewer.pal(3, "Set1")[3]),
                     labels = c(expression(delta[1]), 
                                expression(delta[2]),
                                expression(delta[3]))) +
  labs(x = expression(bar(delta)),
       y = NULL) +
  ggtitle(label = NULL,
          subtitle =  "NPP-SEQ")
plot_hatdelta_onpp_3.1 <- ggplot() +
  geom_density(aes(x = sim3.1$hatdelta[[1]]$hatdelta_onpp, color = "delta1"), linewidth = 0.8) +
  geom_density(aes(x = sim3.1$hatdelta[[2]]$hatdelta_onpp, color = "delta2"), linewidth = 0.8) +
  geom_density(aes(x = sim3.1$hatdelta[[3]]$hatdelta_onpp, color = "delta3"), linewidth = 0.8) +
  scale_color_manual(name = NULL, 
                     values = c("delta1" = RColorBrewer::brewer.pal(3, "Set1")[1], 
                                "delta2" = RColorBrewer::brewer.pal(3, "Set1")[2],
                                "delta3" = RColorBrewer::brewer.pal(3, "Set1")[3]),
                     labels = c(expression(delta[1]), 
                                expression(delta[2]),
                                expression(delta[3]))) +
  labs(x = expression(bar(delta)),
       y = NULL) +
  ggtitle(label = NULL,
          subtitle = "ONPP")
plot_hatdelta_onppseq_3.1 <- ggplot() +
  geom_density(aes(x = sim3.1$hatdelta[[1]]$hatdelta_onppseq, color = "delta1"), linewidth = 0.8) +
  geom_density(aes(x = sim3.1$hatdelta[[2]]$hatdelta_onppseq, color = "delta2"), linewidth = 0.8) +
  geom_density(aes(x = sim3.1$hatdelta[[3]]$hatdelta_onppseq, color = "delta3"), linewidth = 0.8) +
  scale_color_manual(name = NULL, 
                     values = c("delta1" = RColorBrewer::brewer.pal(3, "Set1")[1], 
                                "delta2" = RColorBrewer::brewer.pal(3, "Set1")[2],
                                "delta3" = RColorBrewer::brewer.pal(3, "Set1")[3]),
                     labels = c(expression(delta[1]), 
                                expression(delta[2]),
                                expression(delta[3]))) +
  labs(x = expression(bar(delta)),
       y = NULL) +
  ggtitle(label = NULL,
          subtitle = "ONPP-SEQ")
combined_hatmodels_3.1 <- (plot_hatdelta_npp_3.1 + plot_hatdelta_nppseq_3.1 + plot_hatdelta_onpp_3.1 + plot_hatdelta_onppseq_3.1) + 
  plot_layout(ncol = 4, guides = "collect")
combined_hatmodels_3.1

plot_hatdelta_npp_3.3 <- ggplot() +
  geom_density(aes(x = sim3.3$hatdelta[[1]]$hatdelta_npp, color = "delta1"), linewidth = 0.8) +
  geom_density(aes(x = sim3.3$hatdelta[[2]]$hatdelta_npp, color = "delta2"), linewidth = 0.8) +
  geom_density(aes(x = sim3.3$hatdelta[[3]]$hatdelta_npp, color = "delta3"), linewidth = 0.8) +
  scale_color_manual(name = NULL, 
                     values = c("delta1" = RColorBrewer::brewer.pal(3, "Set1")[1], 
                                "delta2" = RColorBrewer::brewer.pal(3, "Set1")[2],
                                "delta3" = RColorBrewer::brewer.pal(3, "Set1")[3]),
                     labels = c(expression(delta[1]), 
                                expression(delta[2]),
                                expression(delta[3]))) +
  labs(x = expression(bar(delta)),
       y = NULL) +
  ggtitle(label = "Scenario III.III",
          subtitle = "NPP")
plot_hatdelta_nppseq_3.3 <- ggplot() +
  geom_density(aes(x = sim3.3$hatdelta[[1]]$hatdelta_nppseq, color = "delta1"), linewidth = 0.8) +
  geom_density(aes(x = sim3.3$hatdelta[[2]]$hatdelta_nppseq, color = "delta2"), linewidth = 0.8) +
  geom_density(aes(x = sim3.3$hatdelta[[3]]$hatdelta_nppseq, color = "delta3"), linewidth = 0.8) +
  scale_color_manual(name = NULL, 
                     values = c("delta1" = RColorBrewer::brewer.pal(3, "Set1")[1], 
                                "delta2" = RColorBrewer::brewer.pal(3, "Set1")[2],
                                "delta3" = RColorBrewer::brewer.pal(3, "Set1")[3]),
                     labels = c(expression(delta[1]), 
                                expression(delta[2]),
                                expression(delta[3]))) +
  labs(x = expression(bar(delta)),
       y = NULL) +
  ggtitle(label = NULL,
          subtitle =  "NPP-SEQ")
plot_hatdelta_onpp_3.3 <- ggplot() +
  geom_density(aes(x = sim3.3$hatdelta[[1]]$hatdelta_onpp, color = "delta1"), linewidth = 0.8) +
  geom_density(aes(x = sim3.3$hatdelta[[2]]$hatdelta_onpp, color = "delta2"), linewidth = 0.8) +
  geom_density(aes(x = sim3.3$hatdelta[[3]]$hatdelta_onpp, color = "delta3"), linewidth = 0.8) +
  scale_color_manual(name = NULL, 
                     values = c("delta1" = RColorBrewer::brewer.pal(3, "Set1")[1], 
                                "delta2" = RColorBrewer::brewer.pal(3, "Set1")[2],
                                "delta3" = RColorBrewer::brewer.pal(3, "Set1")[3]),
                     labels = c(expression(delta[1]), 
                                expression(delta[2]),
                                expression(delta[3]))) +
  labs(x = expression(bar(delta)),
       y = NULL) +
  ggtitle(label = NULL,
          subtitle = "ONPP")
plot_hatdelta_onppseq_3.3 <- ggplot() +
  geom_density(aes(x = sim3.3$hatdelta[[1]]$hatdelta_onppseq, color = "delta1"), linewidth = 0.8) +
  geom_density(aes(x = sim3.3$hatdelta[[2]]$hatdelta_onppseq, color = "delta2"), linewidth = 0.8) +
  geom_density(aes(x = sim3.3$hatdelta[[3]]$hatdelta_onppseq, color = "delta3"), linewidth = 0.8) +
  scale_color_manual(name = NULL, 
                     values = c("delta1" = RColorBrewer::brewer.pal(3, "Set1")[1], 
                                "delta2" = RColorBrewer::brewer.pal(3, "Set1")[2],
                                "delta3" = RColorBrewer::brewer.pal(3, "Set1")[3]),
                     labels = c(expression(delta[1]), 
                                expression(delta[2]),
                                expression(delta[3]))) +
  labs(x = expression(bar(delta)),
       y = NULL) +
  ggtitle(label = NULL,
          subtitle = "ONPP-SEQ")
combined_hatmodels_3.3 <- (plot_hatdelta_npp_3.3 + plot_hatdelta_nppseq_3.3 + plot_hatdelta_onpp_3.3 + plot_hatdelta_onppseq_3.3) + 
  plot_layout(ncol = 4, guides = "collect")
combined_hatmodels_3.3
