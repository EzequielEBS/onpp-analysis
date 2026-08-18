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
library(gt)
library(gtExtras)

# run auxiliary functions
source("code/aux_fun_sim.R")
source("code/sim_exp/normal/get_summary_normal.r")

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
sim5.2 <- qs_read("results/samples/normal/sim5_2_normal_fixed_var.qs2")

# true value
true_value <- 5

# combine simulations
sim_sces <- list(
  sim1.1, sim1.2, sim1.3,
  sim2.1, sim2.2, sim2.3,
  sim3.1, sim3.2, sim3.3,
  sim4.1, sim4.2, sim4.3,
  sim5.1, sim5.2
)

summary_normal <- get_summary_normal(sim_sces)

# BCI/coverage lists used by the "plot BCI" section below -- these come out
# of get_summary_normal()'s return value; they used to be referenced (as
# bare `bcis_95`, `coverages_95`, etc.) without ever being extracted here,
# which meant this script could only run at all if a same-named object
# happened to be left over in the environment from a previous session (same
# bug as analysis_bin.r had before its fix).
bcis_95 <- summary_normal$bcis_95
bcis_90 <- summary_normal$bcis_90
bcis_80 <- summary_normal$bcis_80
bcis_50 <- summary_normal$bcis_50
coverages_95 <- summary_normal$coverages_95
coverages_90 <- summary_normal$coverages_90
coverages_80 <- summary_normal$coverages_80
coverages_50 <- summary_normal$coverages_50

df_mse_theta <- data.frame(
  # Regime = c(rep("n[0] > n", 12), rep("n[0] < n", 12)),
  Theta0 = c(
    "(4.6, 4.8, 5.0)",
    "(3.6, 3.8, 4.0)",
    "(2.6, 2.8, 3.0)",
    "(4.6, 5.0, 4.8)",
    "(3.6, 4.0, 3.8)",
    "(2.6, 3.0, 2.8)",
    "(5.0, 4.8, 4.6)",
    "(4.0, 3.8, 3.6)",
    "(3.0, 2.8, 2.6)",
    "(5.0, 5.0, 5.0)",
    "(4.0, 4.0, 4.0)",
    "(3.0, 3.0, 3.0)",
    "(3.0, 7.0, 5.0)",
    "(5.0, 7.0, 3.0)"
  ),
  Order       = c(rep("Compatible",    3),
                  rep("Semi-compat.",  3),
                  rep("Incompatible",  3),
                  rep("Neutral",       3),
                  rep("Cancellation effect", 2)
                  ),
  Discrepancy = c(rep(c("None","Small","Large"), 4), rep("Large", 2)),
  NPP         = c(
    lapply(summary_normal$mses, function(x) x[1]) %>% unlist()
  ),
  NPP_SEQ     = c(
    lapply(summary_normal$mses, function(x) x[2]) %>% unlist()
  ),
  ONPP        = c(
    lapply(summary_normal$mses, function(x) x[3]) %>% unlist()
  ),
  ONPP_SEQ    = c(
    lapply(summary_normal$mses, function(x) x[4]) %>% unlist()
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

# apply row-wise coloring anchored at each row's own min/max
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

gtsave(table_mse_theta, "results/summary_tables/normal/mse_theta_table_normal.png", zoom = 3)

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
    lapply(summary_normal$wis_list, function(x) x$hat_wis_npp[1]) %>% unlist()
  ),
  NPP_SEQ     = c(
    lapply(summary_normal$wis_list, function(x) x$hat_wis_nppseq[1]) %>% unlist()
  ),
  ONPP        = c(
    lapply(summary_normal$wis_list, function(x) x$hat_wis_onpp[1]) %>% unlist()
  ),
  ONPP_SEQ    = c(
    lapply(summary_normal$wis_list, function(x) x$hat_wis_onppseq[1]) %>% unlist()
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
gtsave(table_wis_theta, "results/summary_tables/normal/wis_theta_table_normal.png", zoom = 3)

table_mse_wis_theta <- gt_two_column_layout(list(table_mse_theta, table_wis_theta))

htmltools::save_html(table_mse_wis_theta, "table_mse_wis_theta.html")
webshot2::webshot("table_mse_wis_theta.html",
                  "results/summary_tables/normal/mse_wis_theta_table.png",
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
# parametrized by which list it's actually plotting, removes the trap. Same
# fix as analysis_bin.r's build_bci_plot_grid().
#
# The normal model's true theta (5) isn't a probability, so plot_bci() is
# called with model = "Normal" (was "Bernoulli", a leftover from copying
# analysis_bin.r's template -- cosmetic only, since aux_fun_sim.R's
# plot_bci() never actually uses `model` in the plot) and y_range = NULL:
# aux_fun_sim.R's plot_bci() used to hardcode a [0, 1] y-axis window
# (correct for the binomial theta, meaningless here), so every point would
# have been clipped out of a normal-model BCI plot. y_range = NULL lets
# ggplot auto-scale to the actual draws instead.
build_bci_plot_grid <- function(bci_list, coverage_list, prob) {
  plots <- lapply(seq_along(bci_list), function(i) {
    sce <- as.roman(ceiling(i/3))
    sce <- ifelse(i %% 3 == 1, paste0(sce,".I"), ifelse(i %% 3 == 2, paste0(sce,".II"), paste0(sce,".III")))
    title_left <- ggplot() +
      annotate("text", x = 0, y = 0,
               label = sce,
               angle = 90, size = 8) +
      theme_void()
    plot <- plot_bci(bci_list[[i]], sim_sces[[i]], true_value, prob, coverage_list[[i]], sce, "Normal", y_range = NULL)
    plot <- (title_left | plot) + plot_layout(widths = c(0.05, 1))
    if (i == length(bci_list)) {
      plot <- plot & theme(legend.position = "bottom")
    }
    return(plot)
  })
  Reduce(`/`, plots)
}

bci_specs <- list(
  list(bcis = bcis_95, coverages = coverages_95, prob = 0.95, file = "results/figures/normal/bci/bci_95_normal_fixed_var.pdf"),
  list(bcis = bcis_90, coverages = coverages_90, prob = 0.90, file = "results/figures/normal/bci/bci_90_normal_fixed_var.pdf"),
  list(bcis = bcis_80, coverages = coverages_80, prob = 0.80, file = "results/figures/normal/bci/bci_80_normal_fixed_var.pdf"),
  list(bcis = bcis_50, coverages = coverages_50, prob = 0.50, file = "results/figures/normal/bci/bci_50_normal_fixed_var.pdf")
)

for (spec in bci_specs) {
  wrap_plots_bci <- build_bci_plot_grid(spec$bcis, spec$coverages, spec$prob)
  ggsave(spec$file, wrap_plots_bci, width = 21, height = 29, dpi = 300, limitsize = F)
}
