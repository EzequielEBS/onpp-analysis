# =============================================================================
# Logit-normal auxiliary variables (mu, tau, rho, z, e[1:K]) compared across
# the full Order x Discrepancy scenario grid, linear-model (lm) family,
# "varying beta, n0 > n" section -- how does each parameter's posterior
# estimate shift with (i) Discrepancy level (None/Small/Large, pooled
# across Order) and (ii) Order compatibility (Compatible/Semi-compat./
# Incompatible/Neutral, pooled across Discrepancy)?
#
# Direct lm-family port of the bin family's logitnorm_aux_scenario_plot.R --
# the cross-scenario companion to logitnorm_aux_lm_plot.R (which looks at
# ONE scenario's posteriors in detail). Here the unit of comparison is one
# point per (scenario, replicate) -- that replicate's posterior MEDIAN for
# each parameter -- grouped by Discrepancy or Order and shown as boxplots.
#
# IMPORTANT: this script does NOT fit anything itself. prior_comparison_lm_
# sample.R already fits the hierarchical logit-normal prior across every
# scenario in the grid and caches mu/tau/rho/z/e for each (scenario,
# replicate) as aux_logit -- this script just reads that cache and
# re-groups it. Run prior_comparison_lm_sample.R first if the cache doesn't
# exist yet.
# =============================================================================

library(dplyr)
library(tidyverse)
library(qs2)

# Reuses .sce3_theme() and the project's Okabe-Ito color conventions.
source("code/aux_fun_sim.R")

#-------------------------------------------------------------------------------
# Load the cache prior_comparison_lm_sample.R writes.
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
# One row per (scenario, replicate, parameter): that replicate's posterior
# MEDIAN for mu/tau/rho/z/e[1:K], tagged with its scenario's Order and
# Discrepancy. `parameter` is stored as the plotmath string used for facet
# labels (mu, tau, rho, z, epsilon[1..K]) so facet_wrap(labeller =
# label_parsed) can render it directly.
#-------------------------------------------------------------------------------

param_cols <- c("mu", "tau", "rho", "z", sprintf("e[%d]", 1:K))
param_display <- c(mu = "mu", tau = "tau", rho = "rho", z = "z",
                    setNames(sprintf("epsilon[%d]", 1:K), sprintf("e[%d]", 1:K)))

aux_medians <- do.call(rbind, lapply(names(sce_list), function(key) {
  reps <- replicates_by_scenario[[key]]
  sce  <- sce_list[[key]]
  medians <- do.call(rbind, lapply(reps, function(rep) {
    apply(rep$aux_logit[, param_cols, drop = FALSE], 2, median)
  }))
  colnames(medians) <- param_cols
  df <- as.data.frame(medians)
  df$order       <- sce$order
  df$discrepancy <- sce$discrepancy
  df$scenario    <- key
  df
}))

aux_medians_long <- aux_medians %>%
  pivot_longer(cols = all_of(param_cols), names_to = "parameter_key", values_to = "value") %>%
  mutate(parameter   = factor(param_display[parameter_key], levels = unname(param_display)),
         order       = factor(order, levels = order_levels),
         discrepancy = factor(discrepancy, levels = discrepancy_levels))

n_scenarios <- length(unique(aux_medians_long$scenario))

#-------------------------------------------------------------------------------
# (i) By Discrepancy level, pooled across Order. Sequential shading (light
# to dark) since Discrepancy is ordinal (None < Small < Large).
#-------------------------------------------------------------------------------

discrepancy_colors <- c("None" = "#D1E5F0", "Small" = "#4393C3", "Large" = "#053061")

p_aux_by_discrepancy <- ggplot(aux_medians_long, aes(x = discrepancy, y = value, fill = discrepancy)) +
  geom_boxplot(outlier.alpha = 0.3, linewidth = 0.5) +
  facet_wrap(~ parameter, scales = "free_y", labeller = label_parsed, nrow = 2) +
  scale_fill_manual(name = "Discrepancy", values = discrepancy_colors) +
  labs(x = NULL, y = NULL
      #  title = "Logit-normal auxiliary variables by Discrepancy level (lm family)",
      #  subtitle = paste0("Posterior median per replicate, pooled across Order (",
      #                     n_scenarios, " scenarios, ", n_rep_used, " replicates each)")
                        ) +
  .sce3_theme() +
  theme(legend.position = "bottom")

print(p_aux_by_discrepancy)
ggsave("results/figures/lm/logitnorm_aux_lm_by_discrepancy.pdf",
       p_aux_by_discrepancy, width = 12, height = 7)

#-------------------------------------------------------------------------------
# (ii) By Order compatibility, pooled across Discrepancy. Qualitative
# colors (Order is nominal, not ordinal) -- Incompatible picked out in
# vermillion since it's the scenario class that breaks ONPP.
#-------------------------------------------------------------------------------

order_colors <- c("Compatible" = "#0072B2", "Semi-compat." = "#56B4E9",
                   "Incompatible" = "#D55E00", "Neutral" = "#999999")

p_aux_by_order <- ggplot(aux_medians_long, aes(x = order, y = value, fill = order)) +
  geom_boxplot(outlier.alpha = 0.3, linewidth = 0.5) +
  facet_wrap(~ parameter, scales = "free_y", labeller = label_parsed, nrow = 2) +
  scale_fill_manual(name = "Order", values = order_colors) +
  labs(x = NULL, y = NULL
      #  title = "Logit-normal auxiliary variables by Order compatibility (lm family)",
      #  subtitle = paste0("Posterior median per replicate, pooled across Discrepancy (",
      #                     n_scenarios, " scenarios, ", n_rep_used, " replicates each)")
                        ) +
  .sce3_theme() +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 20, hjust = 1))

print(p_aux_by_order)
ggsave("results/figures/lm/logitnorm_aux_lm_by_order.pdf",
       p_aux_by_order, width = 12, height = 7)
