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
source("code/sim_exp/lm/get_summary_lm.r")
load("../../simulated-data-regression/data/true_params.RData")

#-------------------------------------------------------------------------------
# Shared helpers
#-------------------------------------------------------------------------------
#
# The "varying beta" and "varying sigma" sections below each load 12
# scenarios (order I-IV x congruence I-III) for n0 >= n and n0 <= n, then
# build an MSE table and a WIS table styled identically apart from title,
# palette, and whether the Eps/Order/Discrepancy columns are included. That
# used to be ~450 lines of copy-pasted qs_read() calls and gt-table styling
# pipelines; both patterns are shared with simulations_lm.r/analysis_seq_lm.r
# so factoring them out here also makes it obvious at a glance that all four
# summary tables use the same styling rules.

# Loads the 12 sim<order>_<cong>_lm_<suffix>.qs2 files (matching the order
# simulations_lm.r writes them in: order I-IV outer, congruence I-III
# inner) into one list, in the same order the original hand-written
# `list(sim1.1, sim1.2, sim1.3, sim2.1, ...)` blocks used.
load_lm_sims <- function(suffix) {
  sims <- list()
  for (order_i in 1:4) {
    for (cong_j in 1:3) {
      file <- sprintf("results/samples/lm/sim%d_%d_lm_%s.qs2", order_i, cong_j, suffix)
      sims[[length(sims) + 1]] <- qs_read(file)
    }
  }
  sims
}

mse_extractor <- function(summary, model_idx) {
  unlist(lapply(summary$mses, function(x) x[model_idx]))
}
wis_extractor <- function(summary, model_idx) {
  field <- c("hat_wis_npp", "hat_wis_nppseq", "hat_wis_onpp", "hat_wis_onppseq")[model_idx]
  unlist(lapply(summary$wis_list, function(x) mean(x[[field]])))
}

# Builds the 24-row (12 x n0>=n + 12 x n0<=n) metric data frame shared by
# all four summary tables. `extractor(summary, model_idx)` pulls one
# model's values out of a get_summary_lm() result (mse_extractor/
# wis_extractor above); `extra_cols`, if given, is cbind()-ed in before the
# NPP/NPP_SEQ/ONPP/ONPP_SEQ columns (used for the Eps/Order/Discrepancy
# columns in the "varying beta" tables; the "varying sigma" tables omit
# them, as in the original).
build_metric_df <- function(summary_ge, summary_le, extractor, extra_cols = NULL) {
  df <- data.frame(Regime = c(rep("n[0] > n", 12), rep("n[0] < n", 12)))
  if (!is.null(extra_cols)) df <- cbind(df, extra_cols)
  df$NPP      <- c(extractor(summary_ge, 1), extractor(summary_le, 1))
  df$NPP_SEQ  <- c(extractor(summary_ge, 2), extractor(summary_le, 2))
  df$ONPP     <- c(extractor(summary_ge, 3), extractor(summary_le, 3))
  df$ONPP_SEQ <- c(extractor(summary_ge, 4), extractor(summary_le, 4))
  df$Regime <- ifelse(df$Regime == "n[0] > n", "n\u2080 > n", "n\u2080 < n")
  df
}

# Builds and styles one gt summary table (mse_beta/wis_beta/mse_sg/wis_sg
# all use this). `row_palette` is the viridisLite palette function used for
# the per-row (min-to-max) coloring loop -- mako for MSE tables, plasma for
# WIS tables, matching the original.
build_summary_gt_table <- function(df, title, row_palette, include_eps_discrepancy) {
  method_cols <- c("NPP", "NPP_SEQ", "ONPP", "ONPP_SEQ")

  tbl <- df |> gt(groupname_col = "Regime")
  tbl <- if (include_eps_discrepancy) {
    tbl |> cols_label(Eps = "\u03B5", Discrepancy = "Discrepancy",
                       NPP = "NPP", NPP_SEQ = "NPP-SEQ", ONPP = "ONPP", ONPP_SEQ = "ONPP-SEQ")
  } else {
    tbl |> cols_label(NPP = "NPP", NPP_SEQ = "NPP-SEQ", ONPP = "ONPP", ONPP_SEQ = "ONPP-SEQ")
  }

  tbl <- tbl |>
    gt_color_rows(
      columns   = c(NPP, NPP_SEQ, ONPP, ONPP_SEQ),
      palette   = "ggsci::green_material",
      direction = -1
    ) |>
    fmt_number(
      columns  = c(NPP, NPP_SEQ, ONPP, ONPP_SEQ),
      decimals = 3
    ) |>
    tab_header(title = title) |>
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
  for (i in seq_len(nrow(df))) {
    row_vals <- as.numeric(df[i, method_cols])
    tbl <- tbl |>
      data_color(
        columns = all_of(method_cols),
        rows    = i,
        palette = row_palette(100, direction = -1)[30:70],
        domain  = c(min(row_vals), max(row_vals))
      )
  }
  tbl
}

eps_order_discrepancy <- data.frame(
  Eps = rep(c(
    "(1.0, 0.5, 0.0)",
    "(1.5, 1.0, 0.5)",
    "(2.0, 1.5, 1.0)",
    "(1.0, 0.0, 0.5)",
    "(1.5, 0.5, 1.0)",
    "(2.5, 1.5, 2.0)",
    "(0.0, 0.5, 1.0)",
    "(0.5, 1.0, 1.5)",
    "(1.0, 1.5, 2.0)",
    "(0.0, 0.0, 0.0)",
    "(1.0, 1.0, 1.0)",
    "(2.0, 2.0, 2.0)"
  ), 2),
  Order       = rep(c(rep("Compatible",    3),
                      rep("Semi-compat.",  3),
                      rep("Incompatible",  3),
                      rep("Neutral",       3)
                    ), 2),
  Discrepancy = rep(c("None","Small","Large"), 8)
)

#------------------------------------------------------------------------------
# varying beta
#------------------------------------------------------------------------------

betastar <- beta[[1]]
sgstar <- 1

# small size (Order I only: high vs no congruence). This used to be
# assigned to `summary_n0_ge_n`, the same name the full 12-scenario summary
# below gets assigned to -- so it was computed (and its xtable printed) only
# to be immediately overwritten and never read. Given its own name here so
# it's not silently discarded; still unused elsewhere downstream beyond its
# printed xtable, but at least no longer a wasted, invisible computation.
sim1.1_small <- qs_read("results/samples/lm/sim1_1_lm_varying_beta_n0_ge_n_small_size.qs2")
sim1.3_small <- qs_read("results/samples/lm/sim1_3_lm_varying_beta_n0_ge_n_small_size.qs2")
summary_beta_small_size <- get_summary_lm(list(sim1.1_small, sim1.3_small), true_value = c(betastar, sgstar))

sim_sces_n0_ge_n <- load_lm_sims("varying_beta_n0_ge_n")
sim_sces_n0_le_n <- load_lm_sims("varying_beta_n0_le_n")

summary_n0_ge_n <- get_summary_lm(sim_sces_n0_ge_n, true_value = c(betastar, sgstar))
summary_n0_le_n <- get_summary_lm(sim_sces_n0_le_n, true_value = c(betastar, sgstar))

df_mse_beta <- build_metric_df(summary_n0_ge_n, summary_n0_le_n, mse_extractor, eps_order_discrepancy)
table_mse_beta <- build_summary_gt_table(df_mse_beta, "Varying \u03B2", viridisLite::mako, include_eps_discrepancy = TRUE)
gtsave(table_mse_beta, "results/summary_tables/lm/mse_beta_table.png", zoom = 3)

df_wis_beta <- build_metric_df(summary_n0_ge_n, summary_n0_le_n, wis_extractor, eps_order_discrepancy)
table_wis_beta <- build_summary_gt_table(df_wis_beta, "Varying \u03B2", viridisLite::plasma, include_eps_discrepancy = TRUE)
gtsave(table_wis_beta, "results/summary_tables/lm/wis_beta_table.png", zoom = 3)

#------------------------------------------------------------------------------
# varying sigma
#------------------------------------------------------------------------------

sim_sces_n0_ge_n <- load_lm_sims("varying_sigma_n0_ge_n")
sim_sces_n0_le_n <- load_lm_sims("varying_sigma_n0_le_n")

summary_n0_ge_n <- get_summary_lm(sim_sces_n0_ge_n, true_value = c(betastar, sgstar))
summary_n0_le_n <- get_summary_lm(sim_sces_n0_le_n, true_value = c(betastar, sgstar))

df_mse_sg <- build_metric_df(summary_n0_ge_n, summary_n0_le_n, mse_extractor)
table_mse_sg <- build_summary_gt_table(df_mse_sg, "Varying \u03C3", viridisLite::mako, include_eps_discrepancy = FALSE)
gtsave(table_mse_sg, "results/summary_tables/lm/mse_sg_table.png", zoom = 3)

df_wis_sg <- build_metric_df(summary_n0_ge_n, summary_n0_le_n, wis_extractor)
table_wis_sg <- build_summary_gt_table(df_wis_sg, "Varying \u03C3", viridisLite::plasma, include_eps_discrepancy = FALSE)
gtsave(table_wis_sg, "results/summary_tables/lm/wis_sg_table.png", zoom = 3)

#------------------------------------------------------------------------------
# combined layout
#------------------------------------------------------------------------------

table_mse_theta <- gt_two_column_layout(list(table_mse_beta, table_mse_sg))
table_wis_theta <- gt_two_column_layout(list(table_wis_beta, table_wis_sg))

htmltools::save_html(table_mse_theta, "table_mse_theta.html")
htmltools::save_html(table_wis_theta, "table_wis_theta.html")
webshot2::webshot("table_mse_theta.html",
                  "results/summary_tables/lm/mse_theta_table.png",
                  zoom    = 3,
                  vwidth  = 640,   # reduce until white space disappears
                  vheight = 1200)
webshot2::webshot("table_wis_theta.html",
                  "results/summary_tables/lm/wis_theta_table.png",
                  zoom    = 3,
                  vwidth  = 640,   # reduce until white space disappears
                  vheight = 1200)
