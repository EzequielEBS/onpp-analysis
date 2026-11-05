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
# palette, and whether the epsilon/Order/Discrepancy columns are included. That
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

# get_summary_lm() reports the beta/sigma error as an MSE (see
# get_summary_lm.r): for beta (3 coefficients) it's the per-component mean
# squared error, averaged over the 3 coefficients; for sigma (a scalar)
# it's just the mean squared error of that one component, the
# 1-dimensional analogue. WIS is likewise split into WIS_beta (mean of the
# per-coordinate WIS over beta1/beta2/beta3) and WIS_sigma (sigma alone),
# instead of one all-4-parameters-combined figure.
mse_beta_extractor <- function(summary, model_idx) {
  unlist(lapply(summary$mses, function(x) x$mse_beta[model_idx]))
}
mse_sigma_extractor <- function(summary, model_idx) {
  unlist(lapply(summary$mses, function(x) x$mse_sigma[model_idx]))
}
wis_beta_extractor <- function(summary, model_idx) {
  field <- c("hat_wis_npp", "hat_wis_nppseq", "hat_wis_onpp", "hat_wis_onppseq")[model_idx]
  unlist(lapply(summary$wis_list, function(x) mean(x[[field]][1:3])))
}
wis_sigma_extractor <- function(summary, model_idx) {
  field <- c("hat_wis_npp", "hat_wis_nppseq", "hat_wis_onpp", "hat_wis_onppseq")[model_idx]
  unlist(lapply(summary$wis_list, function(x) x[[field]][4]))
}

# Builds the 24-row (12 x n0>=n + 12 x n0<=n) metric data frame shared by
# all four summary tables. For each of the four methods (NPP/NPP-SEQ/
# ONPP/ONPP-SEQ) this now produces two columns -- "<method>_beta" and
# "<method>_sigma" -- via `extractor_beta`/`extractor_sigma`
# (mse_beta_extractor/mse_sigma_extractor or wis_beta_extractor/
# wis_sigma_extractor above). `extra_cols`, if given, is cbind()-ed in
# first (used for the epsilon/Order/Discrepancy columns in the "varying beta"
# tables; the "varying sigma" tables omit them, as in the original).
build_metric_df <- function(summary_ge, summary_le, extractor_beta, extractor_sigma, extra_cols = NULL) {
  df <- data.frame(Regime = c(rep("n[0] > n", 12), rep("n[0] < n", 12)))
  if (!is.null(extra_cols)) df <- cbind(df, extra_cols)
  methods <- c("NPP", "NPP_SEQ", "ONPP", "ONPP_SEQ")
  # columns grouped by parameter (all four methods' beta columns, then all
  # four methods' sigma columns) rather than by method, so the two loops
  # below run separately instead of interleaving beta/sigma per method.
  for (m in seq_along(methods)) {
    df[[paste0(methods[m], "_beta")]] <- c(extractor_beta(summary_ge, m), extractor_beta(summary_le, m))
  }
  for (m in seq_along(methods)) {
    df[[paste0(methods[m], "_sigma")]] <- c(extractor_sigma(summary_ge, m), extractor_sigma(summary_le, m))
  }
  df$Regime <- ifelse(df$Regime == "n[0] > n", "n\u2080 > n", "n\u2080 < n")
  df
}

# Builds and styles one gt summary table (mse_beta/wis_beta/mse_sg/wis_sg
# all use this). `highlight_palette` is the viridisLite palette function
# used for the per-row (min-to-max) gradient -- mako for MSE tables,
# plasma for WIS tables, matching the original; the row's best (smallest)
# cell is additionally marked bold white on top of that gradient. Columns are grouped by parameter
# rather than by method: one spanner for \u03B2 (beta) holding all four
# methods' beta columns, and one spanner for \u03C3 (sigma) holding all
# four methods' sigma columns -- so the two parameters can be compared
# side by side across methods, rather than each method's own beta/sigma
# pair sitting together.
build_summary_gt_table <- function(df, title, highlight_palette, include_eps_discrepancy) {
  methods       <- c("NPP", "NPP_SEQ", "ONPP", "ONPP_SEQ")
  method_labels <- c("NPP", "NPP-SEQ", "ONPP", "ONPP-SEQ")
  beta_cols     <- paste0(methods, "_beta")
  sigma_cols    <- paste0(methods, "_sigma")
  # grouped by parameter: all four methods' beta columns, then all four
  # methods' sigma columns -- NPP_beta, NPP_SEQ_beta, ONPP_beta,
  # ONPP_SEQ_beta, NPP_sigma, NPP_SEQ_sigma, ONPP_sigma, ONPP_SEQ_sigma.
  method_cols <- c(beta_cols, sigma_cols)

  # under each parameter's spanner, the sub-columns are labeled by method
  col_labels <- setNames(as.list(rep(method_labels, 2)), method_cols)
  # NOTE: gt::cols_label() does not reliably apply named `...` labels when
  # `.list` is also supplied in the same call -- only the `.list` entries
  # took effect (that's how the epsilon/Discrepancy headers ended up stuck
  # showing the raw column names while the method-labeled columns, which
  # went through `.list`, rendered correctly). Routing every label through
  # a single `.list` avoids the mixed dots+.list call entirely.
  if (include_eps_discrepancy) {
    col_labels <- c(list(epsilon = "\u03B5", Discrepancy = "Discrepancy"), col_labels)
  }

  tbl <- df |> gt(groupname_col = "Regime") |> cols_label(.list = col_labels)

  tbl <- tbl |>
    fmt_number(
      columns  = all_of(method_cols),
      decimals = 3
    ) |>
    tab_header(title = title) |>
    tab_style(
      style     = cell_text(weight = "bold"),
      locations = cells_row_groups()
    ) |>
    # Rotate the NPP/NPP-SEQ/ONPP/ONPP-SEQ column labels 90 degrees so each
    # one takes up the width of a single narrow column instead of forcing
    # the column wide enough to fit "NPP-SEQ" horizontally -- this is what
    # was making the beta/sigma table too wide to fit the page once each
    # method got its own column (rather than sharing one column split only
    # by beta/sigma). Reads bottom-to-top, the usual convention for
    # rotated table headers (e.g. Excel's default rotation direction).
    tab_style(
      style = css(
        `writing-mode` = "vertical-rl",
        transform      = "rotate(180deg)",
        `white-space`  = "nowrap"
      ),
      locations = cells_column_labels(columns = all_of(method_cols))
    ) |>
    tab_options(
      row_group.font.weight    = "bold",
      row_group.background.color = "#f0f0f0",
      table.font.size          = 11,
      column_labels.padding    = px(8)
    )

  # one spanner for beta (over all four methods) and one for sigma
  tbl <- tbl |>
    tab_spanner(label = "\u03B2", columns = all_of(beta_cols)) |>
    tab_spanner(label = "\u03C3", columns = all_of(sigma_cols))

  # Highlight only the best (smallest) method per row, instead of a
  # continuous per-row color gradient across all eight cells -- see
  # highlight_best_per_row() in aux_fun_sim.R for why. Beta and sigma are
  # on different scales, so the best method is picked *separately* within
  # the beta columns and within the sigma columns (groups =), rather than
  # comparing all eight columns' raw values against each other the way the
  # old single min/max-across-everything domain implicitly did.
  tbl <- highlight_best_per_row(tbl, df, columns = method_cols,
                                groups = list(beta_cols, sigma_cols),
                                palette = highlight_palette)
  tbl
}

eps_order_discrepancy <- data.frame(
  epsilon = rep(c(
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



# NOTE: `sim_sces_*`/`summary_*` below carry a `_beta_`/`_sigma_` tag (unlike
# the earlier shared-name version of this section) specifically so this
# section's variables never collide with the "varying sigma" section's.
# Before this, both sections reused the same `sim_sces_n0_ge_n`/
# `sim_sces_n0_le_n`/`summary_n0_ge_n`/`summary_n0_le_n` names -- fine if
# you always re-source the whole file top to bottom, but re-running just
# one section on its own (e.g. to regenerate one table after a styling
# tweak) silently built that table from whichever section's data happened
# to be sitting in those shared variables, which is exactly what produced
# the "varying beta table shows the same numbers as varying sigma" bug.
# Distinct names make each section runnable independently of the other.
sim_sces_beta_n0_ge_n <- load_lm_sims("varying_beta_n0_ge_n")
sim_sces_beta_n0_le_n <- load_lm_sims("varying_beta_n0_le_n")

summary_beta_n0_ge_n <- get_summary_lm(sim_sces_beta_n0_ge_n, true_value = c(betastar, sgstar))
summary_beta_n0_le_n <- get_summary_lm(sim_sces_beta_n0_le_n, true_value = c(betastar, sgstar))

df_mse_beta <- build_metric_df(summary_beta_n0_ge_n, summary_beta_n0_le_n, mse_beta_extractor, mse_sigma_extractor, eps_order_discrepancy)
table_mse_beta <- build_summary_gt_table(df_mse_beta, "Average MSE \u2014 Varying \u03B2", viridisLite::mako, include_eps_discrepancy = TRUE)
gtsave(table_mse_beta, "results/summary_tables/lm/mse_beta_table.png", zoom = 3)

df_wis_beta <- build_metric_df(summary_beta_n0_ge_n, summary_beta_n0_le_n, wis_beta_extractor, wis_sigma_extractor, eps_order_discrepancy)
table_wis_beta <- build_summary_gt_table(df_wis_beta, "Average WIS \u2014 Varying \u03B2", viridisLite::plasma, include_eps_discrepancy = TRUE)
gtsave(table_wis_beta, "results/summary_tables/lm/wis_beta_table.png", zoom = 3)

#------------------------------------------------------------------------------
# varying sigma
#------------------------------------------------------------------------------

sim_sces_sigma_n0_ge_n <- load_lm_sims("varying_sigma_n0_ge_n")
sim_sces_sigma_n0_le_n <- load_lm_sims("varying_sigma_n0_le_n")

summary_sigma_n0_ge_n <- get_summary_lm(sim_sces_sigma_n0_ge_n, true_value = c(betastar, sgstar))
summary_sigma_n0_le_n <- get_summary_lm(sim_sces_sigma_n0_le_n, true_value = c(betastar, sgstar))

df_mse_sg <- build_metric_df(summary_sigma_n0_ge_n, summary_sigma_n0_le_n, mse_beta_extractor, mse_sigma_extractor)
table_mse_sg <- build_summary_gt_table(df_mse_sg, "Average MSE \u2014 Varying \u03C3", viridisLite::mako, include_eps_discrepancy = FALSE)
gtsave(table_mse_sg, "results/summary_tables/lm/mse_sg_table.png", zoom = 3)

df_wis_sg <- build_metric_df(summary_sigma_n0_ge_n, summary_sigma_n0_le_n, wis_beta_extractor, wis_sigma_extractor)
table_wis_sg <- build_summary_gt_table(df_wis_sg, "Average WIS \u2014 Varying \u03C3", viridisLite::plasma, include_eps_discrepancy = FALSE)
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
                  zoom    = 320/96,
                  vwidth  = 920,   # reduce until white space disappears
                  vheight = 2000)
webshot2::webshot("table_wis_theta.html",
                  "results/summary_tables/lm/wis_theta_table.png",
                  zoom    = 320/96,
                  vwidth  = 920,   # reduce until white space disappears
                  vheight = 2000)
