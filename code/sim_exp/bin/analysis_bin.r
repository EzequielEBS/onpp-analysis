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
library(gt)
library(gtExtras)

# run auxiliary functions
source("code/aux_fun_sim.R")
source("code/sim_exp/bin/get_summary_bin.r")
# scenario definitions (n0/n/a/b/al/bl/alpha/theta0 per scenario) -- shared
# with simulations_bin.R; used below to recover e.g. scenarios[["3_1_all_flat"]]$a
# for the closed-form theta prior plot, instead of relying on untracked globals.
source("code/sim_exp/bin/scenarios_bin.R")

# load sample
sim1.1 <- qs2::qs_read("results/samples/bin/sim1_1_bin.qs2")
sim1.2 <- qs2::qs_read("results/samples/bin/sim1_2_bin.qs2")
sim1.3 <- qs2::qs_read("results/samples/bin/sim1_3_bin.qs2")
sim2.1 <- qs2::qs_read("results/samples/bin/sim2_1_bin.qs2")
sim2.2 <- qs2::qs_read("results/samples/bin/sim2_2_bin.qs2")
sim2.3 <- qs2::qs_read("results/samples/bin/sim2_3_bin.qs2")
sim3.1 <- qs2::qs_read("results/samples/bin/sim3_1_bin.qs2")
sim3.2 <- qs2::qs_read("results/samples/bin/sim3_2_bin.qs2")
sim3.3 <- qs2::qs_read("results/samples/bin/sim3_3_bin.qs2")
sim4.1 <- qs2::qs_read("results/samples/bin/sim4_1_bin.qs2")
sim4.2 <- qs2::qs_read("results/samples/bin/sim4_2_bin.qs2")
sim4.3 <- qs2::qs_read("results/samples/bin/sim4_3_bin.qs2")
sim5.1 <- qs2::qs_read("results/samples/bin/sim5_1_bin.qs2")
sim5.2 <- qs2::qs_read("results/samples/bin/sim5_2_bin.qs2")

sim3.1_eta_flat <- qs2::qs_read("results/samples/bin/sim3_1_bin_eta_flat.qs2")
sim3.2_eta_flat <- qs2::qs_read("results/samples/bin/sim3_2_bin_eta_flat.qs2")
sim3.3_eta_flat <- qs2::qs_read("results/samples/bin/sim3_3_bin_eta_flat.qs2")

sim3.1_eta_gamma_flat <- qs2::qs_read("results/samples/bin/sim3_1_bin_eta_gamma_flat.qs2")
sim3.2_eta_gamma_flat <- qs2::qs_read("results/samples/bin/sim3_2_bin_eta_gamma_flat.qs2")
sim3.3_eta_gamma_flat <- qs2::qs_read("results/samples/bin/sim3_3_bin_eta_gamma_flat.qs2")

sim3.1_all_flat <- qs2::qs_read("results/samples/bin/sim3_1_bin_all_flat.qs2")
sim3.2_all_flat <- qs2::qs_read("results/samples/bin/sim3_2_bin_all_flat.qs2")
sim3.3_all_flat <- qs2::qs_read("results/samples/bin/sim3_3_bin_all_flat.qs2")

samples_prior_sce3.1_all_flat <- qs2::qs_read("results/samples/bin/samples_prior_sce3_1_all_flat.qs2")
samples_prior_sce3.2_all_flat <- qs2::qs_read("results/samples/bin/samples_prior_sce3_2_all_flat.qs2")
samples_prior_sce3.3_all_flat <- qs2::qs_read("results/samples/bin/samples_prior_sce3_3_all_flat.qs2")

# true value
true_value <- 0.5

# combine simulations
sim_sces <- list(
  sim1.1, sim1.2, sim1.3,
  sim2.1, sim2.2, sim2.3,
  sim3.1, sim3.2, sim3.3,
  sim4.1, sim4.2, sim4.3,
  sim5.1, sim5.2
)

sim_sce3 <- list(
  sim3.1, sim3.2, sim3.3,
  sim3.1_eta_flat, sim3.2_eta_flat, sim3.3_eta_flat,
  sim3.1_eta_gamma_flat, sim3.2_eta_gamma_flat, sim3.3_eta_gamma_flat,
  sim3.1_all_flat, sim3.2_all_flat, sim3.3_all_flat
)

summary_bin <- get_summary_bin(sim_sces)
summary_bin_sce3 <- get_summary_bin(sim_sce3)

# BCI/coverage lists used by the "plot BCI" section below -- these come out
# of get_summary_bin()'s return value; they used to be referenced (as bare
# `bcis_95`, `coverages_95`, etc.) without ever being extracted here, which
# meant this script could only run at all if a same-named object happened to
# be left over in the environment from a previous session.
bcis_95 <- summary_bin$bcis_95
bcis_90 <- summary_bin$bcis_90
bcis_80 <- summary_bin$bcis_80
bcis_50 <- summary_bin$bcis_50
coverages_95 <- summary_bin$coverages_95
coverages_90 <- summary_bin$coverages_90
coverages_80 <- summary_bin$coverages_80
coverages_50 <- summary_bin$coverages_50

df_mse_theta <- data.frame(
  # Regime = c(rep("n[0] > n", 12), rep("n[0] < n", 12)),
  Theta0 = c(
    "(0.30, 0.40, 0.50)",
    "(0.20, 0.30, 0.40)",
    "(0.05, 0.15, 0.25)",
    "(0.30, 0.50, 0.40)",
    "(0.20, 0.40, 0.30)",
    "(0.05, 0.25, 0.15)",
    "(0.50, 0.40, 0.30)",
    "(0.40, 0.30, 0.20)",
    "(0.25, 0.15, 0.05)",
    "(0.50, 0.50, 0.50)",
    "(0.30, 0.30, 0.30)",
    "(0.10, 0.10, 0.10)",
    "(0.10, 0.90, 0.50)",
    "(0.50, 0.90, 0.10)"
  ),
  Order       = c(rep("Compatible",    3),
                  rep("Semi-compat.",  3),
                  rep("Incompatible",  3),
                  rep("Neutral",       3),
                  rep("Cancellation effect", 2)
                  ),
  Discrepancy = c(rep(c("None","Small","Large"), 4), rep("Large", 2)),
  NPP         = c(
    lapply(summary_bin$mses, function(x) x[1]) %>% unlist()
  ),
  NPP_SEQ     = c(
    lapply(summary_bin$mses, function(x) x[2]) %>% unlist()
  ),
  ONPP        = c(
    lapply(summary_bin$mses, function(x) x[3]) %>% unlist()
  ),
  ONPP_SEQ    = c(
    lapply(summary_bin$mses, function(x) x[4]) %>% unlist()
  )
)

# df_mse_beta$Regime <- ifelse(df_mse_beta$Regime == "n[0] > n", "n\u2080 > n", "n\u2080 < n")

table_mse_theta <- df_mse_theta |>
  gt(groupname_col = "Order") |>
  cols_label(
    Theta0 = "\u03B8\u2080",
    Discrepancy = "Discrepancy",
    NPP         = "NPP",
    NPP_SEQ     = "NPP-SEQ",
    ONPP        = "ONPP",
    ONPP_SEQ    = "ONPP-SEQ"
  ) |>
  gt_color_rows(
    columns   = c(NPP, NPP_SEQ, ONPP, ONPP_SEQ),
    palette   = "ggsci::green_material",
    direction = -1
  ) |>
  fmt_number(
    columns  = c(NPP, NPP_SEQ, ONPP, ONPP_SEQ),
    decimals = 3
  ) |>
  tab_header(
    # title    = "Average MSE by scenario and method",
    # subtitle = "Varying \u03B2",
    title = "Average MSE"
  ) |>
  tab_spanner(
    label   = "Method",
    columns = c(NPP, NPP_SEQ, ONPP, ONPP_SEQ)
  ) |>
  tab_style(
    style     = cell_text(weight = "bold"),
    locations = cells_row_groups()
  ) |>
  tab_options(
    row_group.font.weight    = "bold",
    row_group.background.color = "#f0f0f0",
    table.font.size          = 11
  )

# apply row-wise coloring anchored at 0
method_cols <- c("NPP", "NPP_SEQ", "ONPP", "ONPP_SEQ")
for (i in seq_len(nrow(df_mse_theta))) {
  row_vals <- as.numeric(df_mse_theta[i, method_cols])
  
  table_mse_theta <- table_mse_theta |>
    data_color(
      columns = all_of(method_cols),
      rows    = i,
      palette = viridisLite::mako(100, direction = -1)[30:70] ,
      domain  = c(min(row_vals), max(row_vals))
    )
}

gtsave(table_mse_theta, "results/summary_tables/bin/mse_theta_table_bin.png", zoom = 3)

df_wis_theta <- data.frame(
  # Theta0 = c(
  #   "(0.30, 0.40, 0.50)",
  #   "(0.20, 0.30, 0.40)",
  #   "(0.05, 0.15, 0.25)",
  #   "(0.30, 0.50, 0.40)",
  #   "(0.20, 0.40, 0.30)",
  #   "(0.05, 0.25, 0.15)",
  #   "(0.50, 0.40, 0.30)",
  #   "(0.40, 0.30, 0.20)",
  #   "(0.25, 0.15, 0.05)",
  #   "(0.50, 0.50, 0.50)",
  #   "(0.30, 0.30, 0.30)",
  #   "(0.10, 0.10, 0.10)",
  #   "(0.10, 0.90, 0.50)"
  # ),
  Order       = c(rep("Compatible",    3),
                  rep("Semi-compat.",  3),
                  rep("Incompatible",  3),
                  rep("Neutral",       3),
                  rep("Cancellation effect", 2)
                  ),
  # Discrepancy = c(rep(c("None","Small","Large"), 4), "Large"),
  NPP         = c(
    lapply(summary_bin$wis_list, function(x) x$hat_wis_npp[1]) %>% unlist()
  ),
  NPP_SEQ     = c(
    lapply(summary_bin$wis_list, function(x) x$hat_wis_nppseq[1]) %>% unlist()
  ),
  ONPP        = c(
    lapply(summary_bin$wis_list, function(x) x$hat_wis_onpp[1]) %>% unlist()
  ),
  ONPP_SEQ    = c(
    lapply(summary_bin$wis_list, function(x) x$hat_wis_onppseq[1]) %>% unlist()
  )
)
table_wis_theta <- df_wis_theta |>
  gt(groupname_col = "Order") |>
  cols_label(
    # Theta0 = "\u03B8\u2080",
    # Discrepancy = "Discrepancy",
    NPP         = "NPP",
    NPP_SEQ     = "NPP-SEQ",
    ONPP        = "ONPP",
    ONPP_SEQ    = "ONPP-SEQ"
  ) |>
  gt_color_rows(
    columns   = c(NPP, NPP_SEQ, ONPP, ONPP_SEQ),
    palette   = "ggsci::green_material",
    direction = -1
  ) |>
  fmt_number(
    columns  = c(NPP, NPP_SEQ, ONPP, ONPP_SEQ),
    decimals = 3
  ) |>
  tab_header(
    title    = "Average WIS"
  ) |>
  tab_spanner(
    label   = "Method",
    columns = c(NPP, NPP_SEQ, ONPP, ONPP_SEQ)
  ) |>
  tab_style(
    style     = cell_text(weight = "bold"),
    locations = cells_row_groups()
  ) |>
  tab_options(
    row_group.font.weight    = "bold",
    row_group.background.color = "#f0f0f0",
    table.font.size          = 11
  )

# apply row-wise coloring anchored at 0
for (i in seq_len(nrow(df_wis_theta))) {
  row_vals <- as.numeric(df_wis_theta[i, method_cols])
  
  table_wis_theta <- table_wis_theta |>
    data_color(
      columns = all_of(method_cols),
      rows    = i,
      palette = viridisLite::plasma(100, direction = -1)[30:70] ,
      domain  = c(min(row_vals), max(row_vals))
    )
}
gtsave(table_wis_theta, "results/summary_tables/bin/wis_theta_table_bin.png", zoom = 3)

table_mse_wis_theta <- gt_two_column_layout(list(table_mse_theta, table_wis_theta))

htmltools::save_html(table_mse_wis_theta, "table_mse_wis_theta.html")
webshot2::webshot("table_mse_wis_theta.html",
                  "results/summary_tables/bin/mse_wis_theta_table.png",
                  zoom    = 3,
                  vwidth  = 590,   # reduce until white space disappears
                  vheight = 1200)

# plot BCI
#
# The four blocks below (95/90/80/50%) used to be near-identical
# copy-pasted lapply()s. Each one's "is this the last row" legend check --
# `if (i == length(bcis_95))` -- was copy-pasted along with the rest of the
# block, so the 90/80/50% versions all compared against `length(bcis_95)`
# instead of their own list. That was harmless only because all four lists
# happen to have the same length (14); factoring the block into one helper,
# parametrized by which list it's actually plotting, removes the trap.
build_bci_plot_grid <- function(bci_list, coverage_list, prob) {
  plots <- lapply(seq_along(bci_list), function(i) {
    sce <- as.roman(ceiling(i/3))
    sce <- ifelse(i %% 3 == 1, paste0(sce,".I"), ifelse(i %% 3 == 2, paste0(sce,".II"), paste0(sce,".III")))
    title_left <- ggplot() +
      annotate("text", x = 0, y = 0,
               label = sce,
               angle = 90, size = 8) +
      theme_void()
    plot <- plot_bci(bci_list[[i]], sim_sces[[i]], true_value, prob, coverage_list[[i]], sce, "Bernoulli")
    plot <- (title_left | plot) + plot_layout(widths = c(0.05, 1))
    if (i == length(bci_list)) {
      plot <- plot & theme(legend.position = "bottom")
    }
    return(plot)
  })
  Reduce(`/`, plots)
}

bci_specs <- list(
  list(bcis = bcis_95, coverages = coverages_95, prob = 0.95, file = "results/figures/bin/bci/bci_95_bin.pdf"),
  list(bcis = bcis_90, coverages = coverages_90, prob = 0.90, file = "results/figures/bin/bci/bci_90_bin.pdf"),
  list(bcis = bcis_80, coverages = coverages_80, prob = 0.80, file = "results/figures/bin/bci/bci_80_bin.pdf"),
  list(bcis = bcis_50, coverages = coverages_50, prob = 0.50, file = "results/figures/bin/bci/bci_50_bin.pdf")
)

for (spec in bci_specs) {
  wrap_plots_bci <- build_bci_plot_grid(spec$bcis, spec$coverages, spec$prob)
  ggsave(spec$file, wrap_plots_bci, width = 21, height = 29, dpi = 300, limitsize = F)
}


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
                          bias_default = unlist(summary_bin_sce3$biass[1:3]),
                          mse_default = unlist(summary_bin_sce3$mses[1:3]),
                          bias_eta_flat = unlist(summary_bin_sce3$biass[4:6]),
                          mse_eta_flat = unlist(summary_bin_sce3$mses[4:6]),
                          bias_eta_gamma_flat = unlist(summary_bin_sce3$biass[7:9]),
                          mse_eta_gamma_flat = unlist(summary_bin_sce3$mses[7:9]),
                          bias_all_flat = unlist(summary_bin_sce3$biass[10:12]),
                          mse_all_flat = unlist(summary_bin_sce3$mses[10:12])
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


# NOTE: this originally tried to compare prior-only draws across all four
# scenario III prior specs (default / eta_flat / eta_gamma_flat / all_flat),
# via `samples_prior_sce3.1`, `samples_prior_sce3.1_eta_flat` and
# `samples_prior_sce3.1_eta_gamma_flat`. Those three were never generated by
# simulations_bin.R (only the "all flat" prior-only draws are ever sampled --
# see samples_prior_sce3.1_all_flat above), so referencing them here always
# errored. Per decision, this now plots only the one prior spec that's
# actually available; extend simulations_bin.R to sample the other three
# scenario III prior specs (prior-only, post = 0) if the full 4-way
# comparison is needed later.
samples_sce3_eta_prior <- list(
    samples_prior_sce3.1_all_flat
)

plots_eta_prior <- plot_eta_prior(samples_sce3_eta_prior)
ggsave("results/figures/bin/prior_eta_priors_bin.pdf",
       plots_eta_prior, width = 10, height = 6, dpi = 300)

samples_sce3_theta_prior <- list(
    rbeta(8000, scenarios[["3_1_all_flat"]]$a, scenarios[["3_1_all_flat"]]$b)
)
plots_theta_prior <- plot_theta_prior(samples_sce3_theta_prior)
plots_theta_prior
ggsave("results/figures/bin/prior_theta_priors_bin.pdf",
       plots_theta_prior, width = 10, height = 3, dpi = 300)


#------------------------------------------------------------------------------
# Cancellation effect scenario
#------------------------------------------------------------------------------

plots_eta_post_canc_eff1_default <- plot_eta_post(sim_sces[[13]], 1)
plots_eta_post_canc_eff1_seq <- plot_eta_post(sim_sces[[13]], 1, seq_norm = T)
plots_eta_post_canc_eff1 <- plots_eta_post_canc_eff1_default / plots_eta_post_canc_eff1_seq
ggsave("results/figures/bin/post_eta_canc_eff1_bin.png", 
       plots_eta_post_canc_eff1, width = 10, height = 6, dpi = 300)
plots_eta_post_canc_eff2_default <- plot_eta_post(sim_sces[[14]], 1)
plots_eta_post_canc_eff2_seq <- plot_eta_post(sim_sces[[14]], 1, seq_norm = T)
plots_eta_post_canc_eff2 <- plots_eta_post_canc_eff2_default / plots_eta_post_canc_eff2_seq
ggsave("results/figures/bin/post_eta_canc_eff2_bin.png",
       plots_eta_post_canc_eff2, width = 10, height = 6, dpi = 300)