# =============================================================================
# Posterior densities of the logistic-normal prior's auxiliary variables
# (mu, tau, rho, z, e[1:K]) -- logitnorm_lm.stan's hierarchical
# compound-symmetry prior on eta, linear-model (lm) family, "varying beta,
# n0 > n" section -- for ONE chosen replicate, plus a summary table of each
# parameter's posterior median, aggregated across replicates.
#
# Direct lm-family port of the bin family's logitnorm_aux_plot.R -- these
# aux variables mean exactly the same thing regardless of which likelihood
# family eta feeds into (Bernoulli for bin, NIG regression for lm), so the
# plots/table below are unchanged in spirit, just reading the lm cache.
#
# IMPORTANT: this script does NOT fit anything itself. It reuses the exact
# same per-(scenario, replicate) fits that prior_comparison_lm_sample.R
# already produces (that script caches mu/tau/rho/z/e alongside its
# eta/theta draws). Run prior_comparison_lm_sample.R first -- once its
# cache file exists (results/samples/lm/prior_comparison_lm_replicates.qs2),
# this script just reads it.
# =============================================================================

library(dplyr)
library(tidyverse)
library(posterior)
library(patchwork)
library(qs2)
library(gt)

# Reuses .sce3_theme() and the project's Okabe-Ito color conventions.
source("code/aux_fun_sim.R")

#-------------------------------------------------------------------------------
# Load the cached fits from prior_comparison_lm_sample.R -- sce_list (for
# Eps/order/discrepancy labels) and replicates_by_scenario (the actual
# per-replicate draws, including aux_logit = mu/tau/rho/z/e).
#-------------------------------------------------------------------------------

cache_file <- "results/samples/lm/prior_comparison_lm_replicates.qs2"
if (!file.exists(cache_file)) {
  stop("Cache not found at ", cache_file, " -- run code/sim_exp/lm/",
       "prior_comparison_lm_sample.R first (it produces this cache as a ",
       "side effect of fitting the scenario grid).")
}
cached <- qs2::qs_read(cache_file)
sce_list <- cached$sce_list
replicates_by_scenario <- cached$replicates_by_scenario

#-------------------------------------------------------------------------------
# Pick one scenario -- "III_I" (Incompatible, Discrepancy None): the same
# cell prior_comparison_lm_sample.R singles out as the one that breaks
# ONPP's order restriction. Swap scenario_key for any other name in
# names(sce_list) to look at a different cell.
#-------------------------------------------------------------------------------

scenario_key <- "III_I"
stopifnot(scenario_key %in% names(sce_list))

sce  <- sce_list[[scenario_key]]
reps <- replicates_by_scenario[[scenario_key]]
K    <- ncol(reps[[1]]$eta_npp)
message(sprintf("Scenario '%s': %d cached replicates available.", scenario_key, length(reps)))

#-------------------------------------------------------------------------------
# One-replicate posterior density plot: mu/tau/rho/z faceted on top (each
# has a different natural scale/support, hence scales = "free_x"), e[1:K]
# overlaid below (same N(0,1)-ish scale, one panel, colored by k).
#
# NOTE on naming: logitnorm_lm.stan's latent shared factor z is a column
# INSIDE `reps[[i]]$aux_logit` (posterior draws) -- unrelated to any data
# field named z elsewhere in the lm pipeline.
#-------------------------------------------------------------------------------

replicate_idx <- 1  # which of the cached replicates to plot -- change freely
stopifnot(replicate_idx >= 1, replicate_idx <= length(reps))

draws_rep <- reps[[replicate_idx]]$aux_logit

hyper_labels <- c(mu = "μ", tau = "τ", rho = "ρ", z = "z")
hyper_long <- do.call(rbind, lapply(names(hyper_labels), function(p) {
  data.frame(value = as.numeric(draws_rep[, p]), parameter = p)
}))
hyper_long$parameter <- factor(hyper_long$parameter, levels = names(hyper_labels))

p_hyper <- ggplot(hyper_long, aes(x = value)) +
  geom_density(fill = "#000000", alpha = 0.15, color = "#000000", linewidth = 0.8) +
  facet_wrap(~ parameter, nrow = 1, scales = "free_x",
             labeller = as_labeller(c(mu = "μ", tau = "τ", rho = "ρ", z = "z"))) +
  labs(x = NULL, y = NULL,
       title = "Hyperparameters: μ, τ, ρ, z") +
  .sce3_theme()

e_cols <- sprintf("e[%d]", 1:K)
e_long <- do.call(rbind, lapply(1:K, function(k) {
  data.frame(value = as.numeric(draws_rep[, e_cols[k]]), component = paste0("e[", k, "]"))
}))
e_long$component <- factor(e_long$component, levels = paste0("e[", 1:K, "]"))

p_e <- ggplot(e_long, aes(x = value, color = component, fill = component)) +
  geom_density(alpha = 0.15, linewidth = 0.8) +
  scale_color_manual(
    name = NULL, breaks = paste0("e[", 1:3, "]"),
    labels = c(expression(epsilon[1]), expression(epsilon[2]), expression(epsilon[3])),
    values = c("e[1]" = "#E69F00", "e[2]" = "#56B4E9", "e[3]" = "#009E73")
  ) +
  scale_fill_manual(
    name = NULL, breaks = paste0("e[", 1:3, "]"),
    labels = c(expression(epsilon[1]), expression(epsilon[2]), expression(epsilon[3])),
    values = c("e[1]" = "#E69F00", "e[2]" = "#56B4E9", "e[3]" = "#009E73")
  ) +
  labs(x = NULL, y = NULL, title = "Idiosyncratic factors: ε[k]") +
  .sce3_theme()

#-------------------------------------------------------------------------------
# eta[k] posterior for the chosen replicate -- the prior's actual draws
# (not the auxiliary mu/tau/rho/z/e[k] that generate them), shown the same
# way as the idiosyncratic factors above for a direct visual comparison.
#-------------------------------------------------------------------------------

eta_rep  <- reps[[replicate_idx]]$eta_logit
eta_cols <- sprintf("eta[%d]", 1:K)
eta_long <- do.call(rbind, lapply(1:K, function(k) {
  data.frame(value = as.numeric(eta_rep[, eta_cols[k]]), component = paste0("η[", k, "]"))
}))
eta_long$component <- factor(eta_long$component, levels = paste0("η[", 1:K, "]"))

p_eta <- ggplot(eta_long, aes(x = value, color = component, fill = component)) +
  geom_density(alpha = 0.15, linewidth = 0.8) +
  scale_color_manual(
    name = NULL, breaks = paste0("η[", 1:3, "]"),
    labels = c(expression(eta[1]), expression(eta[2]), expression(eta[3])),
    values = c("η[1]" = "#E69F00", "η[2]" = "#56B4E9", "η[3]" = "#009E73")
  ) +
  scale_fill_manual(
    name = NULL, breaks = paste0("η[", 1:3, "]"),
    labels = c(expression(eta[1]), expression(eta[2]), expression(eta[3])),
    values = c("η[1]" = "#E69F00", "η[2]" = "#56B4E9", "η[3]" = "#009E73")
  ) +
  labs(x = NULL, y = NULL, title = "Power parameters: η[k]") +
  .sce3_theme()

p_aux <- wrap_plots(p_hyper, p_e, p_eta, ncol = 1, heights = c(1, 1, 1)) +
  plot_annotation(
    # title = bquote("Auxiliary variable posteriors, replicate #" * .(replicate_idx)),
    subtitle = bquote("ε = " * .(sce$eps) *
                         ", Order: " * .(sce$order) * ", Discrepancy: " * .(sce$discrepancy))
  )

print(p_aux)
ggsave(sprintf("results/figures/lm/logitnorm_aux_lm_posterior_%s_rep%d.pdf", scenario_key, replicate_idx),
       p_aux, width = 10, height = 8,
       device = cairo_pdf)

#-------------------------------------------------------------------------------
# Summary table: posterior median of each aux parameter, computed within
# EACH cached replicate, then summarized (median/mean/SD) ACROSS replicates.
#-------------------------------------------------------------------------------

param_cols <- c("mu", "tau", "rho", "z", e_cols)
param_labels <- c("mu" = "μ", "tau" = "τ", "rho" = "ρ", "z" = "z",
                   setNames(paste0("ε[", 1:K, "]"), e_cols))

per_replicate_medians <- do.call(rbind, lapply(reps, function(rep) {
  setNames(as.data.frame(t(apply(rep$aux_logit[, param_cols, drop = FALSE], 2, median))),
           param_cols)
}))

summary_table <- data.frame(
  Parameter = param_labels[param_cols],
  Median = apply(per_replicate_medians, 2, median),
  Mean   = apply(per_replicate_medians, 2, mean),
  SD     = apply(per_replicate_medians, 2, sd)
)
rownames(summary_table) <- NULL
print(summary_table, row.names = FALSE, digits = 3)

table_aux_median <- summary_table |>
  gt() |>
  fmt_number(columns = c(Median, Mean, SD), decimals = 3) |>
  tab_header(
    title = "Posterior median of each auxiliary variable, across replicates",
    subtitle = paste0(scenario_key, ": ε = ", sce$eps, ", ",
                       length(reps), " replicates")
  ) |>
  cols_label(Parameter = "Parameter", Median = "Median", Mean = "Mean", SD = "SD") |>
  tab_options(table.font.size = 11)

gtsave(table_aux_median,
       sprintf("results/summary_tables/lm/logitnorm_aux_lm_median_table_%s.png", scenario_key),
       zoom = 3)

#-------------------------------------------------------------------------------
# eta[k] posterior median across replicates, big figure -- Order doesn't
# matter for the hierarchical logit-normal prior (it has no order
# restriction), so Incompatible and Neutral suffice to see how Discrepancy
# alone shifts the eta posteriors -- one boxplot panel per (Order,
# Discrepancy) cell, faceted in a grid (rows = Order, columns =
# Discrepancy).
#-------------------------------------------------------------------------------

eta_scenario_keys <- c("III_I", "III_II", "III_III", "IV_I", "IV_II", "IV_III")
stopifnot(all(eta_scenario_keys %in% names(sce_list)))

eta_medians_big <- do.call(rbind, lapply(eta_scenario_keys, function(key) {
  reps_k <- replicates_by_scenario[[key]]
  medians <- do.call(rbind, lapply(reps_k, function(rep) {
    apply(rep$eta_logit[, eta_cols, drop = FALSE], 2, median)
  }))
  colnames(medians) <- eta_cols
  df <- as.data.frame(medians)
  df$order       <- sce_list[[key]]$order
  df$discrepancy <- sce_list[[key]]$discrepancy
  df$scenario    <- key
  df
}))

eta_medians_big_long <- eta_medians_big %>%
  pivot_longer(cols = all_of(eta_cols), names_to = "component_key", values_to = "value") %>%
  mutate(component   = factor(paste0("η[", match(component_key, eta_cols), "]"),
                              levels = paste0("η[", 1:K, "]")),
         order       = factor(order, levels = c("Incompatible", "Neutral")),
         discrepancy = factor(discrepancy, levels = c("None", "Small", "Large")))

n_rep_eta <- length(replicates_by_scenario[[eta_scenario_keys[1]]])

p_eta_median_grid <- ggplot(eta_medians_big_long, aes(x = component, y = value, fill = component)) +
  geom_boxplot(outlier.alpha = 0.3, linewidth = 0.5) +
  facet_wrap(order ~ discrepancy, scales = "free_y") +
  scale_fill_manual(
    name = NULL, breaks = paste0("η[", 1:3, "]"),
    labels = c(expression(eta[1]), expression(eta[2]), expression(eta[3])),
    values = c("η[1]" = "#E69F00", "η[2]" = "#56B4E9", "η[3]" = "#009E73")
  ) +
  labs(x = NULL, y = NULL
      #  title = "η[k] posterior median across replicates, by Order x Discrepancy",
      #  subtitle = paste0(n_rep_eta, " replicates per scenario")
      ) +
  .sce3_theme() +
  theme(legend.position = "bottom")

print(p_eta_median_grid)
ggsave("results/figures/lm/logitnorm_eta_lm_median_grid.pdf",
       p_eta_median_grid, width = 12, height = 7,
       device = cairo_pdf)
