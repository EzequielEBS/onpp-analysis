# =============================================================================
# PLOTTING step for the NPP/ONPP/Logit-Normal prior comparison, linear-model
# (lm) family -- "varying beta, n0 > n" scenario section only (12 scenarios:
# Order I-IV x Congruence I-III), matching prior_comparison_lm_sample.R /
# analysis_seq_lm.r's sce_list/df_sce. Builds every figure/table from the
# cache prior_comparison_lm_sample.R writes -- NO Stan sampling happens in
# this script. Run prior_comparison_lm_sample.R first if the cache doesn't
# exist yet.
#
# Compares NPP's default eta prior (i.i.d. Beta, eta_lm.stan) against a
# "couple without ordering" alternative -- the compound-symmetry
# logistic-normal prior (logitnorm_lm.stan) -- and against ONPP
# (order-restricted, gamma_lm.stan) as a reference. Same motivation as the
# bin family's prior_comparison_plot.R (see its header comment) -- ONPP's
# order restriction helps when the assumed order is right, but hurts when
# it's wrong and there's no discrepancy to justify it; the logistic-normal
# priors couple eta[k]'s toward each other without asserting any ranking.
#
# Unlike the bin/normal families, lm's theta is a 4-vector
# (beta0, beta1, beta2, sigma^2) rather than a scalar, so every
# theta-comparison figure here is a 4-panel one (one panel per component),
# and the Order x Discrepancy grid figure is built separately per
# component instead of as a single figure.
#
# Run from the project root.
# =============================================================================

library(dplyr)
library(tidyverse)
library(posterior)
library(patchwork)
library(gt)
library(gtExtras)
library(scoringutils)
library(viridisLite)
library(qs2)

# Reuses .sce3_theme() and the project's Okabe-Ito color conventions; also
# pulls in the ggplot2/dplyr stack above.
source("code/aux_fun_sim.R")

# True (beta0, beta1, beta2, sigma^2) -- same source as analysis_seq_lm.r.
load("../../simulated-data-regression/data/true_params.RData")
true_theta <- c(beta[[1]], 1)
theta_names  <- c("beta0", "beta1", "beta2", "sigma2")
theta_labels <- c(beta0 = expression(beta[0]), beta1 = expression(beta[1]),
                   beta2 = expression(beta[2]), sigma2 = expression(sigma^2))

#-------------------------------------------------------------------------------
# Load the cache prior_comparison_lm_sample.R writes. No fitting happens
# here -- if this errors, run that script first.
#-------------------------------------------------------------------------------

cache_file <- "results/samples/lm/prior_comparison_lm_replicates.qs2"
if (!file.exists(cache_file)) {
  stop("Cache not found at ", cache_file, " -- run code/sim_exp/lm/",
       "prior_comparison_lm_sample.R first (it produces this cache).")
}
cached <- qs2::qs_read(cache_file)
stopifnot(!is.null(cached$sce_list), !is.null(cached$replicates_by_scenario))

sce_list <- cached$sce_list
replicates_by_scenario <- cached$replicates_by_scenario
K <- ncol(replicates_by_scenario[[1]][[1]]$eta_npp)
n_rep_used <- length(replicates_by_scenario[[1]])

order_levels       <- c("Compatible", "Semi-compat.", "Incompatible", "Neutral")
discrepancy_levels <- c("None", "Small", "Large")

#-------------------------------------------------------------------------------
# Plotting helpers -- same method labels/colors/eta_fields as the bin
# family's prior_comparison_plot.R, since they fit the same 3 priors.
#-------------------------------------------------------------------------------

method_labels <- c(npp = "NPP", onpp = "ONPP", logit = "HLNNPP")
method_colors <- c("NPP" = "#0072B2", "ONPP" = "#D55E00", "HLNNPP" = "#E69F00")
eta_fields   <- c(npp = "eta_npp", onpp = "eta_onpp", logit = "eta_logit")
theta_fields <- c(npp = "theta_npp", onpp = "theta_onpp", logit = "theta_logit")
method_order <- c("npp", "onpp", "logit")

make_eta_plot <- function(reps, subtitle, title = NULL) {
  eta_long <- dplyr::bind_rows(lapply(reps, function(rep) {
    dplyr::bind_rows(lapply(names(eta_fields), function(key) {
      mat <- rep[[eta_fields[[key]]]]
      dplyr::bind_rows(lapply(1:K, function(k) {
        data.frame(value = as.numeric(mat[, k]),
                   component = paste0("eta[", k, "]"),
                   method = key)
      }))
    }))
  }))
  eta_long$method <- factor(method_labels[eta_long$method],
                             levels = method_labels[method_order])
  eta_long$component <- factor(eta_long$component, levels = paste0("eta[", 1:K, "]"))

  ggplot(eta_long, aes(x = value, color = component)) +
    geom_density(bounds = c(0, 1), linewidth = 0.8) +
    facet_wrap(~ method, nrow = 1) +
    scale_color_manual(
      name   = NULL,
      breaks = c("eta[1]", "eta[2]", "eta[3]"),
      labels = c(expression(eta[1]), expression(eta[2]), expression(eta[3])),
      values = c("eta[1]" = "#E69F00", "eta[2]" = "#56B4E9", "eta[3]" = "#009E73")
    ) +
    scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
    labs(x = expression(eta), y = NULL, title = title, subtitle = subtitle) +
    .sce3_theme()
}

# One panel per theta component (beta0/beta1/beta2/sigma2), colored by
# method -- the 4-method analogue of plot_theta_lm()'s 2-method
# (Default/SEQ) panels in analysis_seq_lm.r.
make_theta_component_panel <- function(reps, comp_idx, comp_name) {
  comp_long <- dplyr::bind_rows(lapply(reps, function(rep) {
    dplyr::bind_rows(lapply(names(theta_fields), function(key) {
      data.frame(value = as.numeric(rep[[theta_fields[[key]]]][, comp_idx]), method = key)
    }))
  }))
  comp_long$method <- factor(method_labels[comp_long$method],
                              levels = method_labels[method_order])

  ggplot(comp_long, aes(x = value, color = method)) +
    geom_density(linewidth = 0.8) +
    geom_vline(xintercept = true_theta[comp_idx], linetype = "dotted", color = "grey40") +
    scale_color_manual(name = NULL, values = method_colors) +
    labs(x = theta_labels[[theta_names[comp_idx]]], y = NULL) +
    .sce3_theme()
}

make_theta_plot <- function(reps, subtitle, title = NULL) {
  panels <- lapply(1:4, function(j) make_theta_component_panel(reps, j, theta_names[j]))
  wrap_plots(panels, nrow = 1, guides = "collect") +
    plot_annotation(title = title, subtitle = subtitle) &
    theme(legend.position = "bottom")
}

#-------------------------------------------------------------------------------
# Prior-vs-posterior plotting helpers -- facet by stage (Prior/Posterior)
# within each theta component panel, same idea as the bin family's
# make_theta_prior_post_plot() but per-component since theta is 4-D here.
#-------------------------------------------------------------------------------

make_eta_prior_post_plot <- function(prior_draws, post_draws, subtitle, title = NULL) {
  build_long <- function(d, stage) {
    dplyr::bind_rows(lapply(names(eta_fields), function(key) {
      mat <- d[[eta_fields[[key]]]]
      dplyr::bind_rows(lapply(1:K, function(k) {
        data.frame(value = as.numeric(mat[, k]), component = paste0("eta[", k, "]"),
                   method = key, stage = stage)
      }))
    }))
  }
  eta_long <- dplyr::bind_rows(build_long(prior_draws, "Prior"), build_long(post_draws, "Posterior"))
  eta_long$method    <- factor(method_labels[eta_long$method],
                                levels = method_labels[method_order])
  eta_long$component <- factor(eta_long$component, levels = paste0("eta[", 1:K, "]"))
  eta_long$stage     <- factor(eta_long$stage, levels = c("Prior", "Posterior"))

  ggplot(eta_long, aes(x = value, color = component)) +
    geom_density(bounds = c(0, 1), linewidth = 0.8) +
    facet_grid(stage ~ method) +
    scale_color_manual(
      name   = NULL,
      breaks = c("eta[1]", "eta[2]", "eta[3]"),
      labels = c(expression(eta[1]), expression(eta[2]), expression(eta[3])),
      values = c("eta[1]" = "#E69F00", "eta[2]" = "#56B4E9", "eta[3]" = "#009E73")
    ) +
    scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
    labs(x = expression(eta), y = NULL, title = title, subtitle = subtitle) +
    .sce3_theme()
}

make_theta_component_prior_post_panel <- function(prior_draws, post_draws, comp_idx) {
  build_long <- function(d, stage) {
    dplyr::bind_rows(lapply(names(theta_fields), function(key) {
      data.frame(value = as.numeric(d[[theta_fields[[key]]]][, comp_idx]), method = key, stage = stage)
    }))
  }
  comp_long <- dplyr::bind_rows(build_long(prior_draws, "Prior"), build_long(post_draws, "Posterior"))
  comp_long$method <- factor(method_labels[comp_long$method], levels = method_labels[method_order])
  comp_long$stage  <- factor(comp_long$stage, levels = c("Prior", "Posterior"))

  ggplot(comp_long, aes(x = value, color = method)) +
    geom_density(linewidth = 0.8) +
    facet_wrap(~ stage, scales = "free_y") +
    geom_vline(xintercept = true_theta[comp_idx], linetype = "dotted", color = "grey40") +
    scale_color_manual(name = NULL, values = method_colors) +
    labs(x = theta_labels[[theta_names[comp_idx]]], y = NULL) +
    .sce3_theme()
}

make_theta_prior_post_plot <- function(prior_draws, post_draws, subtitle, title = NULL) {
  panels <- lapply(1:4, function(j) make_theta_component_prior_post_panel(prior_draws, post_draws, j))
  wrap_plots(panels, ncol = 2, guides = "collect") +
    plot_annotation(title = title, subtitle = subtitle) &
    theme(legend.position = "bottom")
}

#-------------------------------------------------------------------------------
# Main comparison across the full scenario grid: one figure PER theta
# component (beta0/beta1/beta2/sigma2), each a theta-component posterior
# density faceted Order (rows) x Discrepancy (columns), colored by method --
# the lm analogue of the bin family's single p_theta_grid figure, split into
# 4 since theta isn't a scalar here.
#-------------------------------------------------------------------------------

theta_grid_long <- dplyr::bind_rows(lapply(names(sce_list), function(key) {
  reps <- replicates_by_scenario[[key]]
  sce  <- sce_list[[key]]
  dplyr::bind_rows(lapply(reps, function(rep) {
    dplyr::bind_rows(lapply(names(theta_fields), function(mkey) {
      mat <- rep[[theta_fields[[mkey]]]]
      dplyr::bind_rows(lapply(1:4, function(j) {
        data.frame(value = as.numeric(mat[, j]), component = theta_names[j], method = mkey)
      }))
    }))
  })) %>% mutate(order = sce$order, discrepancy = sce$discrepancy)
}))
theta_grid_long$method <- factor(method_labels[theta_grid_long$method],
                                  levels = method_labels[method_order])
theta_grid_long$order       <- factor(theta_grid_long$order, levels = order_levels)
theta_grid_long$discrepancy <- factor(theta_grid_long$discrepancy, levels = discrepancy_levels)

make_component_grid <- function(comp_name, comp_idx) {
  df <- theta_grid_long %>% filter(component == comp_name)
  ggplot(df, aes(x = value, color = method)) +
    geom_density(linewidth = 0.8) +
    geom_vline(xintercept = true_theta[comp_idx], linetype = "dotted", color = "black", linewidth = 0.5) +
    facet_grid(order ~ discrepancy, scales = "free_x") +
    scale_color_manual(name = NULL, values = method_colors) +
    labs(x = theta_labels[[comp_name]], y = NULL,
         title = paste0(comp_name, " posterior density across the Order x Discrepancy scenario grid"),
         subtitle = paste0("Posterior draws pooled across ", n_rep_used,
                            " replicates per scenario (dotted line = true value)")) +
    .sce3_theme() +
    theme(legend.position = "bottom")
}

p_theta_grid_beta0  <- make_component_grid("beta0", 1)
p_theta_grid_beta1  <- make_component_grid("beta1", 2)
p_theta_grid_beta2  <- make_component_grid("beta2", 3)
p_theta_grid_sigma2 <- make_component_grid("sigma2", 4)

print(p_theta_grid_beta0)
print(p_theta_grid_beta1)
print(p_theta_grid_beta2)
print(p_theta_grid_sigma2)

# ggsave("results/figures/lm/prior_comparison_lm_theta_grid_beta0.pdf",  p_theta_grid_beta0,  width = 12, height = 10)
# ggsave("results/figures/lm/prior_comparison_lm_theta_grid_beta1.pdf",  p_theta_grid_beta1,  width = 12, height = 10)
# ggsave("results/figures/lm/prior_comparison_lm_theta_grid_beta2.pdf",  p_theta_grid_beta2,  width = 12, height = 10)
# ggsave("results/figures/lm/prior_comparison_lm_theta_grid_sigma2.pdf", p_theta_grid_sigma2, width = 12, height = 10)

#-------------------------------------------------------------------------------
# Prior vs. posterior, for the one (scenario, replicate) cell
# prior_comparison_lm_sample.R pre-fit at post = 0 (see that script's
# prior_cell block -- "III_I" / Incompatible, Discrepancy None by default).
# To look at a different cell, change prior_cell_scenario_key/
# prior_cell_rep_idx in prior_comparison_lm_sample.R and re-run it, then
# re-run this script.
#-------------------------------------------------------------------------------

stopifnot(!is.null(cached$prior_cell))
scenario_key <- cached$prior_cell$scenario_key
rep_idx      <- cached$prior_cell$rep_idx
prior_draws  <- cached$prior_cell$draws

sce_chosen <- sce_list[[scenario_key]]
chosen     <- replicates_by_scenario[[scenario_key]][[rep_idx]]
cell_title <- bquote("ε = " * .(sce_chosen$eps) *
                        ", Order: " * .(sce_chosen$order) * ", Discrepancy: " * .(sce_chosen$discrepancy))

p_eta_rep <- make_eta_prior_post_plot(
  prior_draws, chosen, title = cell_title,
  subtitle = paste0("Replicate #", rep_idx, " of ", length(replicates_by_scenario[[scenario_key]]))
)
p_theta_rep <- make_theta_prior_post_plot(
  prior_draws, chosen, title = cell_title,
  subtitle = paste0("Replicate #", rep_idx, " of ", length(replicates_by_scenario[[scenario_key]]),
                     " (dotted line = true value)")
)

print(p_eta_rep)
print(p_theta_rep)

# ggsave(sprintf("results/figures/lm/prior_comparison_lm_eta_%s_rep%d.pdf",   scenario_key, rep_idx), p_eta_rep,   width = 9, height = 7)
# ggsave(sprintf("results/figures/lm/prior_comparison_lm_theta_%s_rep%d.pdf", scenario_key, rep_idx), p_theta_rep, width = 9, height = 9)

#-------------------------------------------------------------------------------
# MSE / WIS summary across the full scenario grid. lm's theta is 4-D, so
# unlike the bin family's scalar table, MSE/WIS here are computed PER theta
# component and then averaged across the 4 components into one number per
# (scenario, method) -- a single-scale summary of "how far off is theta
# overall", at the cost of not showing which component drives a bad score
# (use the per-component grid figures above for that).
#-------------------------------------------------------------------------------

quantile_level <- c(0.025, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.975)

compute_component_mse_wis <- function(reps, method_key, comp_idx) {
  N <- length(reps)
  draws_list <- lapply(reps, function(rep) as.numeric(rep[[theta_fields[[method_key]]]][, comp_idx]))
  hat <- vapply(draws_list, mean, numeric(1))
  mse <- mean((hat - true_theta[comp_idx])^2)
  q <- do.call(rbind, lapply(draws_list, function(d) quantile(d, probs = quantile_level)))
  wis_val <- mean(wis(rep(true_theta[comp_idx], N), q, quantile_level))
  c(mse = mse, wis = wis_val)
}

grid_summary <- dplyr::bind_rows(lapply(names(sce_list), function(key) {
  reps <- replicates_by_scenario[[key]]
  sce  <- sce_list[[key]]

  # Beta (comps 1-3, averaged) and sigma2 (comp 4) kept as SEPARATE columns
  # (not rows) -- one MSE/WIS value per (method, parameter) pair, laid out
  # side by side so each method's β and σ² columns sit next to each other
  # under that method's own column spanner, instead of stacking them as
  # extra rows.
  per_method <- lapply(method_order, function(mkey) {
    comp_scores <- vapply(1:4, function(j) compute_component_mse_wis(reps, mkey, j), numeric(2))
    list(beta = rowMeans(comp_scores[, 1:3, drop = FALSE]), sigma2 = comp_scores[, 4])
  })
  names(per_method) <- method_order

  data.frame(
    Eps = sce$eps, Order = sce$order, Discrepancy = sce$discrepancy,
    MSE_NPP_beta     = per_method$npp$beta[["mse"]],
    MSE_NPP_sigma2   = per_method$npp$sigma2[["mse"]],
    MSE_ONPP_beta    = per_method$onpp$beta[["mse"]],
    MSE_ONPP_sigma2  = per_method$onpp$sigma2[["mse"]],
    MSE_Logit_beta   = per_method$logit$beta[["mse"]],
    MSE_Logit_sigma2 = per_method$logit$sigma2[["mse"]],
    WIS_NPP_beta     = per_method$npp$beta[["wis"]],
    WIS_NPP_sigma2   = per_method$npp$sigma2[["wis"]],
    WIS_ONPP_beta    = per_method$onpp$beta[["wis"]],
    WIS_ONPP_sigma2  = per_method$onpp$sigma2[["wis"]],
    WIS_Logit_beta   = per_method$logit$beta[["wis"]],
    WIS_Logit_sigma2 = per_method$logit$sigma2[["wis"]]
  )
}))
grid_summary$Order       <- factor(grid_summary$Order, levels = order_levels)
grid_summary$Discrepancy <- factor(grid_summary$Discrepancy, levels = discrepancy_levels)
grid_summary <- grid_summary[order(grid_summary$Discrepancy, grid_summary$Order), ]
print(grid_summary, row.names = FALSE, digits = 3)

# ---- scientific-notation helper (copied from prior_comparison_plot.R /
# analysis_bin.r) -- puts each table on a single common power-of-10 scale
# so cells stay near order 1, with the exponent reported once in the title. ----
.superscript_digits <- c(
  "0" = "⁰", "1" = "¹", "2" = "²", "3" = "³", "4" = "⁴",
  "5" = "⁵", "6" = "⁶", "7" = "⁷", "8" = "⁸", "9" = "⁹"
)
to_superscript <- function(n) {
  neg <- n < 0
  digits <- strsplit(as.character(abs(n)), "")[[1]]
  s <- paste(.superscript_digits[digits], collapse = "")
  if (neg) paste0("⁻", s) else s
}
compute_common_scale <- function(x) {
  x <- x[is.finite(x) & x != 0]
  if (length(x) == 0) return(0L)
  round(-floor(log10(median(abs(x)))))
}
add_common_scale <- function(df, cols) {
  exponent <- compute_common_scale(unlist(df[cols]))
  df[cols] <- lapply(df[cols], function(col) col * 10^exponent)
  label <- if (exponent == 0) "" else paste0(" (×10", to_superscript(-exponent), ")")
  list(df = df, exponent = exponent, label = label)
}

beta_cols   <- c("NPP_beta", "ONPP_beta", "Logit_beta")
sigma2_cols <- c("NPP_sigma2", "ONPP_sigma2", "Logit_sigma2")
value_cols  <- c(beta_cols, sigma2_cols)

# ---- MSE table (row-grouped by Discrepancy; β and σ² each get their own
# column spanner, with NPP/ONPP/HLNNPP as the three columns underneath --
# best-cell highlight picks the best method separately within β's columns
# and within σ²'s columns, since the two parameters live on different
# scales, mako gradient) ----
df_mse_theta <- setNames(
  grid_summary[, c("Eps", "Discrepancy", "Order",
                    "MSE_NPP_beta", "MSE_ONPP_beta", "MSE_Logit_beta",
                    "MSE_NPP_sigma2", "MSE_ONPP_sigma2", "MSE_Logit_sigma2")],
  c("Eps", "Discrepancy", "Order", value_cols)
)
mse_scale <- add_common_scale(df_mse_theta, value_cols)
df_mse_theta <- mse_scale$df

table_mse_theta <- df_mse_theta |>
  gt(groupname_col = "Order") |>
  cols_label(Eps = "ε", Discrepancy = "Discrepancy",
             NPP_beta = "NPP", ONPP_beta = "ONPP", Logit_beta = "HLNNPP",
             NPP_sigma2 = "NPP", ONPP_sigma2 = "ONPP", Logit_sigma2 = "HLNNPP") |>
  fmt_number(columns = all_of(value_cols), decimals = 3) |>
  tab_header(title = paste0("Average MSE", mse_scale$label)) |>
  tab_spanner(label = "β", columns = all_of(beta_cols)) |>
  tab_spanner(label = "σ²", columns = all_of(sigma2_cols)) |>
  tab_style(style = cell_text(weight = "bold"), locations = cells_row_groups()) |>
  # Rotate the NPP/ONPP/HLNNPP column labels 90 degrees -- same fix as
  # analysis_lm.R's beta/sigma table, for the same reason: each method now
  # gets its own narrow column (one per β/σ² pair) instead of sharing a
  # wide one, so horizontal headers would force the table too wide.
  tab_style(
    style = css(
      `writing-mode` = "vertical-rl",
      transform      = "rotate(180deg)",
      `white-space`  = "nowrap"
    ),
    locations = cells_column_labels(columns = all_of(value_cols))
  ) |>
  tab_options(row_group.font.weight = "bold", row_group.background.color = "#f0f0f0",
              table.font.size = 11, column_labels.padding = px(8))

table_mse_theta <- highlight_best_per_row(table_mse_theta, df_mse_theta,
                                          columns = value_cols,
                                          groups = list(beta_cols, sigma2_cols),
                                          palette = viridisLite::mako)
gtsave(table_mse_theta, "results/summary_tables/lm/mse_theta_table_prior_comparison_lm.png", zoom = 3)

# ---- WIS table (same column-spanner style, plasma gradient) ----
df_wis_theta <- setNames(
  grid_summary[, c("Order",
                    "WIS_NPP_beta", "WIS_ONPP_beta", "WIS_Logit_beta",
                    "WIS_NPP_sigma2", "WIS_ONPP_sigma2", "WIS_Logit_sigma2")],
  c("Order", value_cols)
)
wis_scale <- add_common_scale(df_wis_theta, value_cols)
df_wis_theta <- wis_scale$df

table_wis_theta <- df_wis_theta |>
  gt(groupname_col = "Order") |>
  cols_label(NPP_beta = "NPP", ONPP_beta = "ONPP", Logit_beta = "HLNNPP",
             NPP_sigma2 = "NPP", ONPP_sigma2 = "ONPP", Logit_sigma2 = "HLNNPP") |>
  fmt_number(columns = all_of(value_cols), decimals = 3) |>
  tab_header(title = paste0("Average WIS", wis_scale$label)) |>
  tab_spanner(label = "β", columns = all_of(beta_cols)) |>
  tab_spanner(label = "σ²", columns = all_of(sigma2_cols)) |>
  tab_style(style = cell_text(weight = "bold"), locations = cells_row_groups()) |>
  tab_style(
    style = css(
      `writing-mode` = "vertical-rl",
      transform      = "rotate(180deg)",
      `white-space`  = "nowrap"
    ),
    locations = cells_column_labels(columns = all_of(value_cols))
  ) |>
  tab_options(row_group.font.weight = "bold", row_group.background.color = "#f0f0f0",
              table.font.size = 11, column_labels.padding = px(8))

table_wis_theta <- highlight_best_per_row(table_wis_theta, df_wis_theta,
                                          columns = value_cols,
                                          groups = list(beta_cols, sigma2_cols),
                                          palette = viridisLite::plasma)
gtsave(table_wis_theta, "results/summary_tables/lm/wis_theta_table_prior_comparison_lm.png", zoom = 3)

# ---- Combined MSE + WIS figure, side by side ----
table_mse_wis_theta <- gt_two_column_layout(list(table_mse_theta, table_wis_theta))
htmltools::save_html(table_mse_wis_theta, "table_mse_wis_prior_comparison_lm.html")
webshot2::webshot("table_mse_wis_prior_comparison_lm.html",
                  "results/summary_tables/lm/mse_wis_prior_comparison_lm_table.png",
                  zoom    = 320/96,
                  vwidth  = 680,
                  vheight = 2000)
