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

sim1.1 <- readRDS("results/samples/lm/sce_I_I.rds")
sim1.2 <- readRDS("results/samples/lm/sce_I_II.rds")
sim1.3 <- readRDS("results/samples/lm/sce_I_III.rds")
sim2.1 <- readRDS("results/samples/lm/sce_II_I.rds")
sim2.2 <- readRDS("results/samples/lm/sce_II_II.rds")
sim2.3 <- readRDS("results/samples/lm/sce_II_III.rds")
sim3.1 <- readRDS("results/samples/lm/sce_III_I.rds")
sim3.2 <- readRDS("results/samples/lm/sce_III_II.rds")
sim3.3 <- readRDS("results/samples/lm/sce_III_III.rds")
sim4.1 <- readRDS("results/samples/lm/sce_IV_I.rds")
sim4.2 <- readRDS("results/samples/lm/sce_IV_II.rds")
sim4.3 <- readRDS("results/samples/lm/sce_IV_III.rds")

sim_sces <- list(
  sim1.1, sim1.2, sim1.3,
  sim2.1, sim2.2, sim2.3,
  sim3.1, sim3.2, sim3.3,
  sim4.1, sim4.2, sim4.3
)

betastar <- c(1, -0.5, 0.5)
sgstar <- 1

hat_par_sces <- lapply(sim_sces, get_hat_par)

plots_sce_lm <- lapply(1:length(hat_par_sces), function(j) {
  plot_sce_lm(j, hat_par_sces)
})
lapply(1:4, function(i) {
  sce_plots <- unlist(lapply(1:3, function(j) {
    plots_sce_lm[[(i-1)*3 + j]]
  }), recursive = FALSE)
  combined_plot <- patchwork::wrap_plots(sce_plots, ncol = 1)
  sce <- as.roman(i)
  ggsave(filename = paste0("results/figures/lm/boxplot_sce_", sce, "_lm.pdf"),
         plot = combined_plot,
         width = 10, height = 13)
})
  

# Compute MSE
true_value <- c(betastar, sgstar)
mses <- lapply(hat_par_sces, function(sim) {
  colSums(
    do.call(rbind, lapply(seq_along(sim$hattheta), function(i) {
      colMeans((sim$hattheta[[i]] - true_value[i])^2)
    }))
  )
})

biass <- lapply(hat_par_sces, function(sim) {
  do.call(rbind, lapply(seq_along(sim$hattheta), function(i) {
    colMeans((sim$hattheta[[i]] - true_value[i]))
  }))
})

quantile_level <- c(0.05, 0.1, 0.2, 0.5, 0.8, 0.9, 0.95)
wis_list <- lapply(sim_sces, function(sim) {
  compute_wis_lm(sim, quantile_level, true_value)
})

# ------------------------------------------------------------------------------
# Summary results
# ------------------------------------------------------------------------------

summary_table <- data.frame(scenario = unlist(lapply(1:length(sim_sces), 
                                                     function(i) {
                                                       sce <- as.roman(ceiling(i/3))
                                                       sce <- ifelse(i %% 3 == 1, paste0(sce,".I"), ifelse(i %% 3 == 2, paste0(sce,".II"), paste0(sce,".III")))
                                                       rep(sce,4)
                                                     })),
                            model = rep(c("NPP", "NPP-SEQ", "ONPP", "ONPP-SEQ"),length(sim_sces)),
                            MSE = unlist(mses),
                            bias_beta1 = unlist(lapply(biass, function(x) x[1,])),
                            bias_beta2 = unlist(lapply(biass, function(x) x[2,])),
                            bias_beta3 = unlist(lapply(biass, function(x) x[3,])),
                            bias_sg = unlist(lapply(biass, function(x) x[4,])),
                            WIS_beta1 = unlist(lapply(wis_list, function(x) c(x$hat_wis_npp[1],
                                                                              x$hat_wis_nppseq[1],
                                                                              x$hat_wis_onpp[1],
                                                                              x$hat_wis_onppseq[1]))),
                            WIS_beta2 = unlist(lapply(wis_list, function(x) c(x$hat_wis_npp[2],
                                                                              x$hat_wis_nppseq[2],
                                                                              x$hat_wis_onpp[2],
                                                                              x$hat_wis_onppseq[2]))),
                            WIS_beta3 = unlist(lapply(wis_list, function(x) c(x$hat_wis_npp[3],
                                                                              x$hat_wis_nppseq[3],
                                                                              x$hat_wis_onpp[3],
                                                                              x$hat_wis_onppseq[3]))),
                            WIS_sg = unlist(lapply(wis_list, function(x) c(x$hat_wis_npp[4],
                                                                            x$hat_wis_nppseq[4],
                                                                            x$hat_wis_onpp[4],
                                                                            x$hat_wis_onppseq[4])))
)

row.names(summary_table) <- NULL

print(xtable::xtable(summary_table %>% 
                       mutate(
                         MSE = formatC(MSE, format = "f", digits = 3),
                         bias_beta1 = formatC(bias_beta1, format = "f", digits = 3),
                         bias_beta2 = formatC(bias_beta2, format = "f", digits = 3),
                         bias_beta3 = formatC(bias_beta3, format = "f", digits = 3),
                         bias_sg = formatC(bias_sg, format = "f", digits = 3),
                         WIS_beta1 = formatC(WIS_beta1, format = "f", digits = 3),
                         WIS_beta2 = formatC(WIS_beta2, format = "f", digits = 3),
                         WIS_beta3 = formatC(WIS_beta3, format = "f", digits = 3),
                         WIS_sg = formatC(WIS_sg, format = "f", digits = 3)
                       )
), 
type = "latex", include.rownames = FALSE)
