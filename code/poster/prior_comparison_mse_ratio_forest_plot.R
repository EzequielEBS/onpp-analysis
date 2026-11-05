# =============================================================================
# MSE-ratio forest plot for the NPP vs. logistic-normal-prior comparison
# (lm family, prior_comparison_lm_plot.R's cache): ONPP and HLNNPP
# (hierarchical logistic-normal "couple without ordering" prior), each
# shown as MSE relative to NPP -- the same style/layout as
# code/poster/combined_mse_ratio_forest_plot.R, but for this one
# comparison instead of across all three experiments, and with 2 panels
# (beta, sigma^2) instead of 4.
#
# Does NOT re-fit anything -- reuses the exact same cache
# prior_comparison_lm_plot.R reads (results/samples/lm/
# prior_comparison_lm_replicates.qs2) and the same per-component MSE
# computation it uses (compute_component_mse_wis()), just reshaped into a
# ratio-to-NPP long data frame instead of a gt table. Run
# prior_comparison_lm_sample.R first if that cache doesn't exist yet.
#
# Layout: a narrow unlabeled strip on the left carries the Order group
# names (Compatible / Semi-compat. / Incompatible / Neutral) ONCE, glued
# via patchwork to a 2-panel facet_wrap (beta, sigma^2). Both panels share
# the same row_id levels, so with the default (fixed, not free_y) y scale,
# the Discrepancy labels (none/small/large) are drawn only on the
# leftmost panel, not repeated on both.
# =============================================================================

library(dplyr)
library(tidyverse)
library(posterior)
library(patchwork)
library(scoringutils)
library(qs2)

source("code/aux_fun_sim.R")

load("../../simulated-data-regression/data/true_params.RData")
true_theta <- c(beta[[1]], 1)

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

order_levels       <- c("Compatible", "Semi-compat.", "Incompatible", "Neutral")
discrepancy_levels <- c("None", "Small", "Large")

theta_fields <- c(npp = "theta_npp", onpp = "theta_onpp", logit = "theta_logit")
method_order <- c("npp", "onpp", "logit")

quantile_level <- c(0.025, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.975)

compute_component_mse <- function(reps, method_key, comp_idx) {
  draws_list <- lapply(reps, function(rep) as.numeric(rep[[theta_fields[[method_key]]]][, comp_idx]))
  hat <- vapply(draws_list, mean, numeric(1))
  mean((hat - true_theta[comp_idx])^2)
}

#-------------------------------------------------------------------------------
# MSE per (scenario, method, component), beta (comps 1-3, averaged) and
# sigma2 (comp 4) kept separate -- same split prior_comparison_lm_plot.R
# uses for its gt table -- then reshaped into ratio-to-NPP.
#-------------------------------------------------------------------------------

grid_summary <- dplyr::bind_rows(lapply(names(sce_list), function(key) {
  reps <- replicates_by_scenario[[key]]
  sce  <- sce_list[[key]]

  per_method <- lapply(method_order, function(mkey) {
    comp_mse <- vapply(1:4, function(j) compute_component_mse(reps, mkey, j), numeric(1))
    list(beta = mean(comp_mse[1:3]), sigma2 = comp_mse[4])
  })
  names(per_method) <- method_order

  data.frame(
    Order = sce$order, Discrepancy = sce$discrepancy,
    MSE_NPP_beta     = per_method$npp$beta,
    MSE_NPP_sigma2   = per_method$npp$sigma2,
    MSE_ONPP_beta    = per_method$onpp$beta,
    MSE_ONPP_sigma2  = per_method$onpp$sigma2,
    MSE_HLN_beta     = per_method$logit$beta,
    MSE_HLN_sigma2   = per_method$logit$sigma2
  )
}))

ratio_df <- grid_summary %>%
  mutate(
    ONPP_beta   = MSE_ONPP_beta   / MSE_NPP_beta,
    HLN_beta    = MSE_HLN_beta    / MSE_NPP_beta,
    ONPP_sigma2 = MSE_ONPP_sigma2 / MSE_NPP_sigma2,
    HLN_sigma2  = MSE_HLN_sigma2  / MSE_NPP_sigma2
  ) %>%
  select(Order, Discrepancy, ONPP_beta, HLN_beta, ONPP_sigma2, HLN_sigma2)

panel_labels <- c(beta = "β", sigma2 = "σ²")

df_beta <- ratio_df %>%
  select(Order, Discrepancy, ONPP = ONPP_beta, `HLN-NPP` = HLN_beta) %>%
  pivot_longer(cols = c("ONPP", "HLN-NPP"), names_to = "method", values_to = "ratio") %>%
  mutate(panel = panel_labels[["beta"]])
df_sigma2 <- ratio_df %>%
  select(Order, Discrepancy, ONPP = ONPP_sigma2, `HLN-NPP` = HLN_sigma2) %>%
  pivot_longer(cols = c("ONPP", "HLN-NPP"), names_to = "method", values_to = "ratio") %>%
  mutate(panel = panel_labels[["sigma2"]])

all_df <- bind_rows(df_beta, df_sigma2)
all_df$panel  <- factor(all_df$panel, levels = unname(panel_labels))
all_df$method <- factor(all_df$method, levels = c("ONPP", "HLN-NPP"))
all_df$Order       <- factor(all_df$Order, levels = order_levels)
all_df$Discrepancy <- factor(all_df$Discrepancy, levels = discrepancy_levels)

row_order <- all_df %>%
  distinct(Order, Discrepancy) %>%
  arrange(Order, Discrepancy) %>%
  mutate(row_id = paste(Order, Discrepancy, sep = "__"))

all_df$row_id <- factor(paste(all_df$Order, all_df$Discrepancy, sep = "__"),
                         levels = rev(row_order$row_id))

method_colors <- c("ONPP" = "#0072B2", "HLN-NPP" = "#009E73")

#-------------------------------------------------------------------------------
# Alternating background bands, one per Order group.
#-------------------------------------------------------------------------------

order_band_df <- row_order %>%
  mutate(y_num = as.numeric(factor(row_id, levels = rev(row_order$row_id)))) %>%
  group_by(Order) %>%
  summarise(ymin = min(y_num) - 0.5, ymax = max(y_num) + 0.5, .groups = "drop") %>%
  mutate(fill_group = as.numeric(factor(Order, levels = order_levels)) %% 2)

# Thin separator line between every consecutive pair of scenario rows
# (not just between Order groups, which the bands above already set off).
scenario_sep_y <- seq(1.5, nrow(row_order) - 0.5, by = 1)

# Bold separator line between consecutive ORDER groups (Compatible /
# Semi-compat. / Incompatible / Neutral / Cancellation) -- order_band_df's
# rows are already in order_levels order (group_by/summarise on a factor
# sorts by level), so each group's own ymax is the boundary below it;
# drop the last one, which is the outer edge rather than a boundary.
order_sep_y <- order_band_df$ymax[-1]

#-------------------------------------------------------------------------------
# Main 2-panel plot.
#-------------------------------------------------------------------------------

p_main <- ggplot(all_df, aes(x = ratio, y = row_id)) +
  geom_rect(data = order_band_df,
            aes(xmin = -Inf, xmax = Inf, ymin = ymin, ymax = ymax, fill = factor(fill_group)),
            inherit.aes = FALSE, alpha = 0.6) +
  scale_fill_manual(values = c("0" = "white", "1" = "grey93"), guide = "none") +
  geom_hline(yintercept = scenario_sep_y, color = "grey80", linewidth = 0.3) +
  geom_hline(yintercept = order_sep_y, color = "black", linewidth = 0.9) +
  geom_vline(xintercept = 1, color = "black", linewidth = 0.7) +
  geom_point(aes(color = method), size = 3.4) +
  scale_color_manual(values = method_colors, name = NULL) +
  scale_x_log10(breaks = c(0.5, 0.67, 1, 1.5, 2), labels = c("0.5", "0.67", "1", "1.5", "2")) +
  scale_y_discrete(labels = function(x) tolower(sub("^.*__", "", x))) +
  facet_wrap(~ panel, nrow = 1) +
  labs(x = "MSE relative to NPP", y = NULL) +
  .sce3_theme() +
  theme(
    legend.position    = "top",
    legend.text        = element_text(size = 12),
    strip.text         = element_text(face = "bold.italic", size = 16, hjust = 0),
    panel.grid.major.y = element_blank(),
    axis.text.y        = element_text(size = 12),
    axis.text.x        = element_text(size = 11),
    axis.title          = element_text(size = 14)
  )

#-------------------------------------------------------------------------------
# Narrow left-hand strip carrying the Order group names ONCE.
#-------------------------------------------------------------------------------

order_label_df <- order_band_df %>% mutate(y_mid = (ymin + ymax) / 2)

p_order_labels <- ggplot(order_label_df, aes(x = 0, y = y_mid, label = Order)) +
  geom_text(hjust = 1, size = 4.6, color = "#2C5F8A", fontface = "bold") +
  coord_cartesian(xlim = c(-1, 0), ylim = c(0.5, nrow(row_order) + 0.5), clip = "off") +
  theme_void() 

p_combined <- p_order_labels + p_main + plot_layout(widths = c(0.14, 1)) +
  theme(
    legend.position = "none"
  )

print(p_combined)
ggsave("results/figures/poster/prior_comparison_mse_ratio_forest_plot.pdf",
       p_combined, width = 10.2, height = 3.5, dev = cairo_pdf)
