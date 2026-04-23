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
load("results/samples/bin/sim5_1_bin.RData")

load("results/samples/bin/sim3_1_bin_eta_flat.RData")
load("results/samples/bin/sim3_2_bin_eta_flat.RData")
load("results/samples/bin/sim3_3_bin_eta_flat.RData")
load("results/samples/bin/sim3_1_bin_eta_gamma_flat.RData")
load("results/samples/bin/sim3_2_bin_eta_gamma_flat.RData")
load("results/samples/bin/sim3_3_bin_eta_gamma_flat.RData")
load("results/samples/bin/sim3_1_bin_all_flat.RData")
load("results/samples/bin/sim3_2_bin_all_flat.RData")
load("results/samples/bin/sim3_3_bin_all_flat.RData")

# combine simulations
sim_sces <- list(
  sim1.1, sim1.2, sim1.3,
  sim2.1, sim2.2, sim2.3,
  sim3.1, sim3.2, sim3.3,
  sim4.1, sim4.2, sim4.3,
  sim5.1
)

sim_sce3 <- list(
  sim3.1, sim3.2, sim3.3,
  sim3.1_eta_flat, sim3.2_eta_flat, sim3.3_eta_flat,
  sim3.1_eta_gamma_flat, sim3.2_eta_gamma_flat, sim3.3_eta_gamma_flat,
  sim3.1_all_flat, sim3.2_all_flat, sim3.3_all_flat
)

# Compute MSE
true_value <- 0.5
mses <- lapply(sim_sces, function(sim) {
  colMeans((sim$hattheta - true_value)^2)
})
mses_sce3 <- lapply(sim_sce3, function(sim) {
  colMeans((sim$hattheta - true_value)^2)
})

biass <- lapply(sim_sces, function(sim) {
  colMeans(sim$hattheta) - true_value
})
biass_sce3 <- lapply(sim_sce3, function(sim) {
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
# Comparison of scenarios 3
# ------------------------------------------------------------------------------

comp_table <- data.frame(scenario = unlist(lapply(7:9, 
                                    function(i) {
                                      sce <- as.roman(ceiling(i/3))
                                      sce <- ifelse(i %% 3 == 1, paste0(sce,".I"), ifelse(i %% 3 == 2, paste0(sce,".II"), paste0(sce,".III")))
                                      rep(sce,4)
                                    })),
                          model = rep(c("NPP", "NPP-SEQ", "ONPP", "ONPP-SEQ"),3),
                          bias_default = unlist(biass_sce3[1:3]),
                          mse_default = unlist(mses_sce3[1:3]),
                          bias_eta_flat = unlist(biass_sce3[4:6]),
                          mse_eta_flat = unlist(mses_sce3[4:6]),
                          bias_eta_gamma_flat = unlist(biass_sce3[7:9]),
                          mse_eta_gamma_flat = unlist(mses_sce3[7:9]),
                          bias_all_flat = unlist(biass_sce3[10:12]),
                          mse_all_flat = unlist(mses_sce3[10:12])
)
colnames(comp_table) <- c("Scenario", 
                          "Prior", 
                          "Bias_Default", 
                          "MSE_Default",
                          "Bias_Eta_Flat",
                          "MSE_Eta_Flat",
                          "Bias_Eta_Gamma_Flat",
                          "MSE_Eta_Gamma_Flat",
                          "Bias_All_Flat",
                          "MSE_All_Flat"
                          )

print(xtable::xtable(comp_table %>% 
                       mutate(
                         Bias_Default = formatC(Bias_Default, format = "f", digits = 3),
                         MSE_Default = formatC(MSE_Default, format = "f", digits = 3),
                         Bias_Eta_Flat = formatC(Bias_Eta_Flat, format = "f", digits = 3),
                         MSE_Eta_Flat = formatC(MSE_Eta_Flat, format = "f", digits = 3),
                         Bias_Eta_Gamma_Flat = formatC(Bias_Eta_Gamma_Flat, format = "f", digits = 3),
                         MSE_Eta_Gamma_Flat = formatC(MSE_Eta_Gamma_Flat, format = "f", digits = 3),
                         Bias_All_Flat = formatC(Bias_All_Flat, format = "f", digits = 3),
                         MSE_All_Flat = formatC(MSE_All_Flat, format = "f", digits = 3)
                       )),
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
# Plot distributions for a given simulation
# ------------------------------------------------------------------------------

sim <- 5  # choose which simulation to plot

sce3_priors <- list(
  list(
    sim3.1, sim3.1_eta_flat, sim3.1_eta_gamma_flat, sim3.1_all_flat
  ),
  list(
    sim3.2, sim3.2_eta_flat, sim3.2_eta_gamma_flat, sim3.2_all_flat
  ),
  list(
     sim3.3, sim3.3_eta_flat, sim3.3_eta_gamma_flat, sim3.3_all_flat
  )
)

plots_sce3_eta <- lapply(1:3, function(i) {
  plot_sce3_eta(i, sim = sim, sim_sces = sce3_priors)
})
ggsave("results/figures/bin/post_eta_sce_III_priors_bin.pdf", 
        wrap_plots(plots_sce3_eta, ncol = 1), width = 10, height = 12, dpi = 300)

plots_sce3_theta <- lapply(1:3, function(i) {
  plot_sce3_theta(i, sim = sim, sim_sces = sce3_priors)
})
ggsave("results/figures/bin/post_theta_sce_III_priors_bin.pdf", 
       wrap_plots(plots_sce3_theta, ncol = 1), width = 10, height = 8, dpi = 300)


samples_sce3_eta_prior <- list(
    samples_prior_sce3.1,
    samples_prior_sce3.1_eta_flat,
    samples_prior_sce3.1_eta_gamma_flat,
    samples_prior_sce3.1_all_flat
)

plots_eta_prior <- plot_eta_prior(samples_sce3_eta_prior)
ggsave("results/figures/bin/prior_eta_priors_bin.pdf", 
       plots_eta_prior, width = 10, height = 6, dpi = 300)

samples_sce3_theta_prior <- list(
    rbeta(8000, sce3.1$a, sce3.1$b),
    rbeta(8000, sce3.1_eta_flat$a, sce3.1_eta_flat$b),
    rbeta(8000, sce3.1_eta_gamma_flat$a, sce3.1_eta_gamma_flat$b),
    rbeta(8000, sce3.1_all_flat$a, sce3.1_all_flat$b)
)
plots_theta_prior <- plot_theta_prior(samples_sce3_theta_prior)
plots_theta_prior
ggsave("results/figures/bin/prior_theta_priors_bin.pdf", 
       plots_theta_prior, width = 10, height = 3, dpi = 300)

# plots_sce_eta <- lapply(1:length(sim_sces), function(i) {
#   plot_sim_delta(i, sim = sim, sim_sces = sim_sces)
# })

# lapply(1:4, function(i) {
#   sce_plots <- lapply(1:3, function(j) {
#     plots_sce_eta[[(i-1)*3 + j]]
#   })
#   combined_plot <- Reduce(`/`, sce_plots)
#   sce <- as.roman(i)
#   ggsave(paste0("results/figures/bin/post_eta_sce_", sce, "_bin.pdf"), combined_plot, width = 10, height = 12, dpi = 300)
# })

# plots_sce_theta <- lapply(1:length(sim_sces), function(i) {
#   plot_sim_theta(i, sim = sim, sim_sces = sim_sces)
# })

# lapply(1:4, function(i) {
#   sce_plots <- lapply(1:3, function(j) {
#     plots_sce_theta[[(i-1)*3 + j]]
#   })
#   combined_plot <- Reduce(`/`, sce_plots)
#   sce <- as.roman(i)
#   ggsave(paste0("results/figures/bin/post_theta_sce_", sce, "_bin.pdf"), combined_plot, width = 10, height = 12, dpi = 300)
# })


# plots_sce_eta_os <- lapply(1:length(sim_sces_os), function(i) {
#   plot_sim_delta(i, sim = sim, sim_sces = sim_sces_os)
# })
# lapply(3, function(i) {
#   sce_plots <- lapply(1:3, function(j) {
#     plots_sce_eta_os[[(i-1)*3 + j]]
#   })
#   combined_plot <- Reduce(`/`, sce_plots)
#   sce <- as.roman(i)
#   ggsave(paste0("results/figures/bin/post_eta_sce_", sce, "_bin_os.pdf"), combined_plot, width = 10, height = 12, dpi = 300)
# })
