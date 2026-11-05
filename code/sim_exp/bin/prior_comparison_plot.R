# =============================================================================
# PLOTTING step for the NPP/ONPP/Logit-Normal prior comparison across the
# full Order x Discrepancy scenario grid (Bernoulli experiments, Section
# 5.1). Builds every figure/table from the cache prior_comparison_sample.R
# writes -- NO Stan sampling happens in this script. Run
# prior_comparison_sample.R first if the cache doesn't exist yet.
#
# Compares NPP's default eta prior (i.i.d. Beta) against a "couple without
# ordering" alternative -- the compound-symmetry logistic-normal prior
# (logitnorm_bin.stan) -- and against ONPP (order-restricted,
# gamma_bin.stan) as a reference.
#
# Motivation (see the manuscript's Section 6 desiderata and the
# Incompatible/None failure mode discussed in Sections 5.1/5.3): ONPP's order
# restriction helps whenever the assumed order is right or genuine
# discrepancy is present, but hurts badly when the order is WRONG and there
# is no discrepancy to justify the extra shrinkage. The logistic-normal
# prior below couples the eta[k]'s toward each other -- reproducing part of
# ONPP's shrinkage behaviour -- WITHOUT asserting any ranking across the K
# historical datasets, so it can't be "wrong" about an order the way ONPP
# can.
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

#-------------------------------------------------------------------------------
# Load the cache prior_comparison_sample.R writes. No fitting happens here
# -- if this errors, run that script first.
#-------------------------------------------------------------------------------

cache_file <- "results/samples/bin/prior_comparison_replicates.qs2"
if (!file.exists(cache_file)) {
  stop("Cache not found at ", cache_file, " -- run code/sim_exp/bin/",
       "prior_comparison_sample.R first (it produces this cache).")
}
cached <- qs2::qs_read(cache_file)
stopifnot(!is.null(cached$sce_list), !is.null(cached$replicates_by_scenario))

sce_list <- cached$sce_list
replicates_by_scenario <- cached$replicates_by_scenario
K <- length(sce_list[[1]]$n0)
true_theta <- sce_list[[1]]$theta  # shared across every scenario in the grid
n_rep_used <- length(replicates_by_scenario[[1]])

order_levels      <- c("Compatible", "Semi-compatible", "Incompatible", "Neutral")
discrepancy_levels <- c("None", "Small", "Large")

#-------------------------------------------------------------------------------
# Plotting helpers, built to take ANY subset of a scenario's replicates -- a
# single one (list(reps[[i]])) or all of them (pooled). Boundary-corrected
# KDE on (0,1), same convention as plot_sce3_eta()/plot_eta_prior_post() in
# aux_fun_sim.R.
#-------------------------------------------------------------------------------

method_labels <- c(npp = "NPP", onpp = "ONPP", logit = "HLNNPP")
method_colors <- c("NPP" = "#0072B2", "ONPP" = "#D55E00", "HLNNPP" = "#E69F00")
eta_fields <- c(npp = "eta_npp", onpp = "eta_onpp", logit = "eta_logit")
method_order <- c("npp", "onpp", "logit")

make_eta_plot <- function(reps, subtitle, title = NULL) {
  eta_long <- do.call(rbind, lapply(reps, function(rep) {
    do.call(rbind, lapply(names(eta_fields), function(key) {
      mat <- rep[[eta_fields[[key]]]]
      do.call(rbind, lapply(1:K, function(k) {
        data.frame(value = as.numeric(mat[, k]),
                   component = paste0("eta[", k, "]"),
                   method = key)
      }))
    }))
  }))
  eta_long$method <- factor(method_labels[eta_long$method],
                             levels = method_labels[method_order])
  eta_long$component <- factor(eta_long$component, levels = paste0("eta[", 1:3, "]"))

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

make_theta_plot <- function(reps, subtitle, title = NULL) {
  theta_long <- do.call(rbind, lapply(reps, function(rep) {
    rbind(
      data.frame(value = rep$theta_npp,   method = "npp"),
      data.frame(value = rep$theta_onpp,  method = "onpp"),
      data.frame(value = rep$theta_logit, method = "logit")
    )
  }))
  theta_long$method <- factor(method_labels[theta_long$method],
                               levels = method_labels[method_order])

  ggplot(theta_long, aes(x = value, color = method)) +
    geom_density(bounds = c(0, 1), linewidth = 0.9) +
    scale_color_manual(name = NULL, values = method_colors) +
    scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
    labs(x = expression(theta), y = NULL, title = title, subtitle = subtitle) +
    .sce3_theme()
}

#-------------------------------------------------------------------------------
# Prior-vs-posterior plotting helpers -- same idea as make_eta_plot()/
# make_theta_plot() above, but each takes a PRIOR draws list and a POSTERIOR
# draws list together and facets by stage (matching the two-panel Prior/
# Posterior layout of theta_prior_post_npp.pdf / eta_prior_post_npp.pdf in
# aux_fun_sim.R -- those compare Default vs. SEQ within each stage; these
# compare all 3 eta-priors within each stage instead).
#-------------------------------------------------------------------------------

make_eta_prior_post_plot <- function(prior_draws, post_draws, subtitle, title = NULL) {
  build_long <- function(d, stage) {
    do.call(rbind, lapply(names(eta_fields), function(key) {
      mat <- d[[eta_fields[[key]]]]
      do.call(rbind, lapply(1:K, function(k) {
        data.frame(value = as.numeric(mat[, k]), component = paste0("eta[", k, "]"),
                   method = key, stage = stage)
      }))
    }))
  }
  eta_long <- rbind(build_long(prior_draws, "Prior"), build_long(post_draws, "Posterior"))
  eta_long$method    <- factor(method_labels[eta_long$method],
                                levels = method_labels[method_order])
  eta_long$component <- factor(eta_long$component, levels = paste0("eta[", 1:3, "]"))
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

make_theta_prior_post_plot <- function(prior_draws, post_draws, subtitle, title = NULL) {
  build_long <- function(d, stage) {
    rbind(
      data.frame(value = d$theta_npp,   method = "npp",   stage = stage),
      data.frame(value = d$theta_onpp,  method = "onpp",  stage = stage),
      data.frame(value = d$theta_logit, method = "logit", stage = stage)
    )
  }
  theta_long <- rbind(build_long(prior_draws, "Prior"), build_long(post_draws, "Posterior"))
  theta_long$method <- factor(method_labels[theta_long$method],
                               levels = method_labels[method_order])
  theta_long$stage  <- factor(theta_long$stage, levels = c("Prior", "Posterior"))

  ggplot(theta_long, aes(x = value, color = method)) +
    geom_density(bounds = c(0, 1), linewidth = 0.9) +
    facet_wrap(~ stage) +
    scale_color_manual(name = NULL, values = method_colors) +
    scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
    labs(x = expression(theta), y = NULL, title = title, subtitle = subtitle) +
    .sce3_theme()
}

#-------------------------------------------------------------------------------
# Main comparison across the full scenario grid: theta posterior density,
# pooled across each scenario's replicates, faceted Order (rows) x
# Discrepancy (columns), colored by method. This is the primary "compare
# across scenarios" figure -- an eta-level version isn't included here since
# faceting by component AND method AND order AND discrepancy at once (a
# 4-way facet) stops being legible; use make_eta_plot()/
# make_eta_prior_post_plot() on one scenario at a time (see below) instead.
#-------------------------------------------------------------------------------

theta_grid_long <- do.call(rbind, lapply(names(sce_list), function(key) {
  reps <- replicates_by_scenario[[key]]
  sce  <- sce_list[[key]]
  do.call(rbind, lapply(reps, function(rep) {
    rbind(
      data.frame(value = rep$theta_npp,   method = "npp"),
      data.frame(value = rep$theta_onpp,  method = "onpp"),
      data.frame(value = rep$theta_logit, method = "logit")
    )
  })) %>% mutate(order = sce$order, discrepancy = sce$discrepancy)
}))
theta_grid_long$method <- factor(method_labels[theta_grid_long$method],
                                  levels = method_labels[method_order])
theta_grid_long$order       <- factor(theta_grid_long$order, levels = order_levels)
theta_grid_long$discrepancy <- factor(theta_grid_long$discrepancy, levels = discrepancy_levels)

p_theta_grid <- ggplot(theta_grid_long, aes(x = value, color = method)) +
  geom_density(bounds = c(0, 1), linewidth = 0.8) +
  geom_vline(xintercept = true_theta, linetype = "dotted", color = "black", linewidth = 0.5) +
  facet_grid(order ~ discrepancy) +
  scale_color_manual(name = NULL, values = method_colors) +
  scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
  labs(x = expression(theta), y = NULL,
       title = "theta posterior density across the Order x Discrepancy scenario grid",
       subtitle = paste0("Posterior draws pooled across ", n_rep_used,
                          " replicates per scenario (dotted line = true θ = ", true_theta, ")")) +
  .sce3_theme()

print(p_theta_grid)

# ggsave("results/figures/bin/prior_comparison_theta_grid.pdf", p_theta_grid, width = 12, height = 10)

#-------------------------------------------------------------------------------
# Prior vs. posterior, for the one (scenario, replicate) cell
# prior_comparison_sample.R pre-fit at post = 0 (see that script's
# prior_cell block -- "Incompatible_None" by default, the scenario that
# breaks ONPP's order restriction). To look at a different cell, change
# prior_cell_scenario_key/prior_cell_rep_idx in prior_comparison_sample.R
# and re-run it, then re-run this script.
#-------------------------------------------------------------------------------

stopifnot(!is.null(cached$prior_cell))
scenario_key <- cached$prior_cell$scenario_key
rep_idx      <- cached$prior_cell$rep_idx
prior_draws  <- cached$prior_cell$draws

sce_chosen <- sce_list[[scenario_key]]
chosen     <- replicates_by_scenario[[scenario_key]][[rep_idx]]
cell_title <- bquote(theta[0] == .(paste0("(", paste(formatC(sce_chosen$theta0, format = "f", digits = 2), collapse = ", "), ")")) *
                        ", Order: " * .(sce_chosen$order) * ", Discrepancy: " * .(sce_chosen$discrepancy))

p_eta_rep <- make_eta_prior_post_plot(
  prior_draws, chosen, title = cell_title,
  subtitle = paste0("Replicate #", rep_idx, " of ", length(replicates_by_scenario[[scenario_key]]))
)
p_theta_rep <- make_theta_prior_post_plot(
  prior_draws, chosen, title = cell_title,
  subtitle = paste0("Replicate #", rep_idx, " of ", length(replicates_by_scenario[[scenario_key]]),
                     " (dotted line = true θ = ", sce_chosen$theta, ")")
)

print(p_eta_rep)
print(p_theta_rep)

# ggsave(sprintf("results/figures/bin/prior_comparison_eta_%s_rep%d.pdf",   scenario_key, rep_idx), p_eta_rep,   width = 9, height = 7)
# ggsave(sprintf("results/figures/bin/prior_comparison_theta_%s_rep%d.pdf", scenario_key, rep_idx), p_theta_rep, width = 7, height = 5)

#-------------------------------------------------------------------------------
# MSE / WIS summary across the full scenario grid (theta posterior mean vs.
# truth, and the quantile-based Weighted Interval Score vs. truth) -- the
# same two quantities Figure 3 / the "Average MSE"/"Average WIS" tables
# report in analysis_bin.r, computed here for the 3-method eta-prior
# comparison (NPP / ONPP / Logit-Normal) across this script's 12-scenario
# Order x Discrepancy grid instead of analysis_bin.r's 14-scenario table.
#
# WIS uses scoringutils::wis() the same way compute_wis() does in
# aux_fun_sim.R: per-replicate theta quantiles at the standard 9 levels,
# scored against the true theta, then averaged across replicates.
#-------------------------------------------------------------------------------

quantile_level <- c(0.025, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.975)

grid_summary <- do.call(rbind, lapply(names(sce_list), function(key) {
  reps <- replicates_by_scenario[[key]]
  sce  <- sce_list[[key]]
  N    <- length(reps)

  hattheta <- do.call(rbind, lapply(reps, function(rep) {
    data.frame(npp = mean(rep$theta_npp), onpp = mean(rep$theta_onpp),
               logit = mean(rep$theta_logit))
  }))
  mse <- colMeans((hattheta - sce$theta)^2)

  q_npp   <- do.call(rbind, lapply(reps, function(rep) quantile(rep$theta_npp,   probs = quantile_level)))
  q_onpp  <- do.call(rbind, lapply(reps, function(rep) quantile(rep$theta_onpp,  probs = quantile_level)))
  q_logit <- do.call(rbind, lapply(reps, function(rep) quantile(rep$theta_logit, probs = quantile_level)))
  wis_npp   <- mean(wis(rep(sce$theta, N), q_npp,   quantile_level))
  wis_onpp  <- mean(wis(rep(sce$theta, N), q_onpp,  quantile_level))
  wis_logit <- mean(wis(rep(sce$theta, N), q_logit, quantile_level))

  data.frame(
    Theta0 = paste0("(", paste(formatC(sce$theta0, format = "f", digits = 2), collapse = ", "), ")"),
    Order = sce$order, Discrepancy = sce$discrepancy,
    MSE_NPP = mse[["npp"]], MSE_ONPP = mse[["onpp"]], MSE_Logit = mse[["logit"]],
    WIS_NPP = wis_npp, WIS_ONPP = wis_onpp, WIS_Logit = wis_logit
  )
}))
grid_summary$Order       <- factor(grid_summary$Order, levels = order_levels)
grid_summary$Discrepancy <- factor(grid_summary$Discrepancy, levels = discrepancy_levels)
grid_summary <- grid_summary[order(grid_summary$Order, grid_summary$Discrepancy), ]
print(grid_summary, row.names = FALSE, digits = 3)

# ---- scientific-notation helper (copied from analysis_bin.r, not in
# aux_fun_sim.R) -- puts each table on a single common power-of-10 scale so
# cells stay near order 1, with the exponent reported once in the title. ----
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

method_cols3 <- c("NPP", "ONPP", "Logit_Normal")

# ---- MSE table (row-grouped by Order, mako gradient, best-cell highlight --
# same visual style as analysis_bin.r's table_mse_theta) ----
df_mse_theta <- setNames(
  grid_summary[, c("Theta0", "Order", "Discrepancy", "MSE_NPP", "MSE_ONPP", "MSE_Logit")],
  c("Theta0", "Order", "Discrepancy", method_cols3)
)
mse_scale <- add_common_scale(df_mse_theta, method_cols3)
df_mse_theta <- mse_scale$df

table_mse_theta <- df_mse_theta |>
  gt(groupname_col = "Order") |>
  cols_label(Theta0 = "θ₀", Discrepancy = "Discrepancy",
             NPP = "NPP", ONPP = "ONPP", Logit_Normal = "HLNNPP") |>
  fmt_number(columns = all_of(method_cols3), decimals = 3) |>
  tab_header(title = paste0("Average MSE", mse_scale$label)) |>
  tab_spanner(label = "Method", columns = all_of(method_cols3)) |>
  tab_style(style = cell_text(weight = "bold"), locations = cells_row_groups()) |>
  tab_options(row_group.font.weight = "bold", row_group.background.color = "#f0f0f0",
              table.font.size = 11)

table_mse_theta <- highlight_best_per_row(table_mse_theta, df_mse_theta,
                                          columns = method_cols3,
                                          palette = viridisLite::mako)
gtsave(table_mse_theta, "results/summary_tables/bin/mse_theta_table_prior_comparison.png", zoom = 3)

# ---- WIS table (same style, plasma gradient -- matching table_wis_theta) ----
df_wis_theta <- setNames(
  grid_summary[, c("Theta0", "Order", "Discrepancy", "WIS_NPP", "WIS_ONPP", "WIS_Logit")],
  c("Theta0", "Order", "Discrepancy", method_cols3)
)
wis_scale <- add_common_scale(df_wis_theta, method_cols3)
df_wis_theta <- wis_scale$df

table_wis_theta <- df_wis_theta |>
  gt(groupname_col = "Order") |>
  cols_label(Theta0 = "θ₀", Discrepancy = "Discrepancy",
             NPP = "NPP", ONPP = "ONPP", Logit_Normal = "HLNNPP") |>
  fmt_number(columns = all_of(method_cols3), decimals = 3) |>
  tab_header(title = paste0("Average WIS", wis_scale$label)) |>
  tab_spanner(label = "Method", columns = all_of(method_cols3)) |>
  tab_style(style = cell_text(weight = "bold"), locations = cells_row_groups()) |>
  tab_options(row_group.font.weight = "bold", row_group.background.color = "#f0f0f0",
              table.font.size = 11)

table_wis_theta <- highlight_best_per_row(table_wis_theta, df_wis_theta,
                                          columns = method_cols3,
                                          palette = viridisLite::plasma)
gtsave(table_wis_theta, "results/summary_tables/bin/wis_theta_table_prior_comparison.png", zoom = 3)

# ---- Combined MSE + WIS figure, side by side (gt_two_column_layout, same
# as analysis_bin.r's table_mse_wis_theta) ----
table_mse_wis_theta <- gt_two_column_layout(list(table_mse_theta, table_wis_theta))
htmltools::save_html(table_mse_wis_theta, "table_mse_wis_prior_comparison.html")
webshot2::webshot("table_mse_wis_prior_comparison.html",
                  "results/summary_tables/bin/mse_wis_prior_comparison_table.png",
                  zoom    = 3,
                  vwidth  = 1000,
                  vheight = 900)
