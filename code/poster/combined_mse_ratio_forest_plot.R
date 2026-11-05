# =============================================================================
# Combined cross-experiment MSE-ratio forest plot: NPP-SEQ / ONPP / ONPP-SEQ
# MSE, each shown relative to NPP's MSE -- NPP is this project's baseline
# "joint" normalization prior (eta_bin.stan / eta_lm.stan: a single shared
# power parameter across all K historical datasets, no per-dataset or
# sequential discounting) -- across all three simulation experiments:
# Bernoulli theta (analysis_bin.r), Normal-fixed-variance theta
# (analysis_normal.r), and linear regression beta/sigma^2 (analysis_lm.R,
# "varying beta, n0 > n" setting only, per the user's scoping -- lm has no
# "varying sigma" or Cancellation-effect analogue here).
#
# Does NOT re-fit anything -- reuses each family's own get_summary_*()
# helper on the exact same cached .qs2 replicate files each analysis_*.R/
# analysis_*.r script already reads, and works off the RAW (unscaled) MSE
# numbers: any common power-of-10 scaling used elsewhere (see
# add_common_scale() in those scripts) cancels out of a ratio, so there's
# nothing to undo here.
#
# Two-part layout: a narrow unlabeled strip on the left carries the Order
# group names (Compatible / Semi-compat. / Incompatible / Neutral /
# Cancellation) ONCE, glued via patchwork to a 4-panel facet_wrap (one
# panel per experiment). All four panels share the same row_id levels, so
# with the default (fixed, not free_y) y scale, ggplot only draws the
# Discrepancy labels (none/small/large, or A/B for Cancellation) on the
# leftmost panel -- not repeated across the other three -- while the
# x-axis (the ratio scale) stays shared/fixed across all four panels too.
# =============================================================================

library(dplyr)
library(tidyverse)
library(qs2)
library(patchwork)
library(scoringutils)

source("code/aux_fun_sim.R")

method_names <- c("NPP", "NPP_SEQ", "ONPP", "ONPP_SEQ")

# Turns a list of per-scenario length-4 MSE vectors (model order NPP /
# NPP-SEQ / ONPP / ONPP-SEQ, as every get_summary_*() helper in this
# project returns -- positional, not name-based, matching how
# analysis_bin.r/analysis_normal.r/analysis_lm.R already extract them)
# into one ratio-to-NPP long data frame tagged with this panel's name and
# the row labels (Order/Discrepancy) passed in.
mse_ratio_long <- function(mses_list, row_labels, panel_name) {
  mse_mat <- do.call(rbind, lapply(mses_list, function(x) as.numeric(x)))
  colnames(mse_mat) <- method_names
  ratio_mat <- mse_mat / mse_mat[, "NPP"]
  df <- as.data.frame(ratio_mat[, c("NPP_SEQ", "ONPP", "ONPP_SEQ")])
  df <- cbind(row_labels, df)
  df$panel <- panel_name
  df %>% pivot_longer(cols = c("NPP_SEQ", "ONPP", "ONPP_SEQ"),
                       names_to = "method", values_to = "ratio")
}

#-------------------------------------------------------------------------------
# Bernoulli (bin): theta MSE -- 14 scenarios (12 Order x Discrepancy + 2
# Cancellation-effect sub-scenarios, labeled A/B here rather than the
# "Large" Discrepancy analysis_bin.r otherwise tags them with, matching
# the reference figure).
#-------------------------------------------------------------------------------

source("code/sim_exp/bin/get_summary_bin.r")

bin_sim_names <- c("sim1_1", "sim1_2", "sim1_3", "sim2_1", "sim2_2", "sim2_3",
                    "sim3_1", "sim3_2", "sim3_3", "sim4_1", "sim4_2", "sim4_3",
                    "sim5_1", "sim5_2")
sim_sces_bin <- lapply(bin_sim_names, function(nm) {
  qs2::qs_read(sprintf("results/samples/bin/%s_bin.qs2", nm))
})
summary_bin <- get_summary_bin(sim_sces_bin, true_value = 0.5)

scenario_row_labels <- data.frame(
  Order       = c(rep("Compatible", 3), rep("Semi-compat.", 3),
                   rep("Incompatible", 3), rep("Neutral", 3),
                   rep("Cancellation", 2)),
  Discrepancy = c(rep(c("none", "small", "large"), 4), "A", "B")
)
df_bin <- mse_ratio_long(summary_bin$mses, scenario_row_labels, "Bernoulli: theta")

#-------------------------------------------------------------------------------
# Normal (fixed-variance): theta MSE -- same 14-scenario Order x
# Discrepancy (+ Cancellation A/B) layout as bin.
#-------------------------------------------------------------------------------

source("code/sim_exp/normal/get_summary_normal.r")

sim_sces_normal <- lapply(bin_sim_names, function(nm) {
  qs2::qs_read(sprintf("results/samples/normal/%s_normal_fixed_var.qs2", nm))
})
summary_normal <- get_summary_normal(sim_sces_normal, true_value = 5)

df_normal <- mse_ratio_long(summary_normal$mses, scenario_row_labels, "Normal: theta")

#-------------------------------------------------------------------------------
# Linear regression: beta and sigma^2 MSE -- "varying beta, n0 > n"
# setting only (per the user's scoping). 12 scenarios (Order x
# Discrepancy), no Cancellation-effect analogue -- lm has none (see
# analysis_lm.R's header comment).
#-------------------------------------------------------------------------------

source("code/sim_exp/lm/get_summary_lm.r")
load("../../simulated-data-regression/data/true_params.RData")
betastar <- beta[[1]]
sgstar <- 1

sim_sces_lm_beta_ge <- list()
for (order_i in 1:4) {
  for (cong_j in 1:3) {
    file <- sprintf("results/samples/lm/sim%d_%d_lm_varying_beta_n0_ge_n.qs2", order_i, cong_j)
    sim_sces_lm_beta_ge[[length(sim_sces_lm_beta_ge) + 1]] <- qs_read(file)
  }
}
summary_lm_beta_ge <- get_summary_lm(sim_sces_lm_beta_ge, true_value = c(betastar, sgstar))

lm_row_labels <- data.frame(
  Order       = c(rep("Compatible", 3), rep("Semi-compat.", 3),
                   rep("Incompatible", 3), rep("Neutral", 3)),
  Discrepancy = rep(c("none", "small", "large"), 4)
)
mses_beta  <- lapply(summary_lm_beta_ge$mses, function(x) x$mse_beta)
mses_sigma <- lapply(summary_lm_beta_ge$mses, function(x) x$mse_sigma)
df_lm_beta  <- mse_ratio_long(mses_beta,  lm_row_labels, "Regression: beta")
df_lm_sigma <- mse_ratio_long(mses_sigma, lm_row_labels, "Regression: sigma^2")

# lm has no Cancellation scenario -- two placeholder rows (ratio = NA, so
# nothing is drawn there) so all four panels still share the same 14-row
# y-axis and line up; "not defined" is labeled in-panel over them below.
lm_cancellation_placeholder <- expand.grid(
  Order = "Cancellation", Discrepancy = c("A", "B"),
  method = c("NPP_SEQ", "ONPP", "ONPP_SEQ"), stringsAsFactors = FALSE
)
lm_cancellation_placeholder$ratio <- NA_real_

df_lm_beta  <- bind_rows(df_lm_beta,  cbind(lm_cancellation_placeholder, panel = "Regression: beta"))
df_lm_sigma <- bind_rows(df_lm_sigma, cbind(lm_cancellation_placeholder, panel = "Regression: sigma^2"))

#-------------------------------------------------------------------------------
# Combine and set up the shared row ordering (bottom-to-top factor levels,
# since ggplot's first discrete-y level renders at the BOTTOM): Cancellation
# at the very bottom, Compatible at the very top, matching the reference
# figure; within each Order group, none/small/large (or A/B) read top-down.
#-------------------------------------------------------------------------------

all_df <- bind_rows(df_bin, df_normal, df_lm_beta, df_lm_sigma)

panel_order  <- c("Bernoulli: theta", "Normal: theta", "Regression: beta", "Regression: sigma^2")
panel_labels <- c("Bernoulli: theta"     = "Bernoulli: θ",
                   "Normal: theta"       = "Normal: θ",
                   "Regression: beta"    = "Regression: β",
                   "Regression: sigma^2" = "Regression: σ²")
all_df$panel <- factor(panel_labels[all_df$panel], levels = unname(panel_labels[panel_order]))

order_levels <- c("Compatible", "Semi-compat.", "Incompatible", "Neutral", "Cancellation")
all_df$Order       <- factor(all_df$Order, levels = order_levels)
all_df$Discrepancy <- factor(all_df$Discrepancy, levels = c("none", "small", "large", "A", "B"))

row_order <- all_df %>%
  distinct(Order, Discrepancy) %>%
  arrange(Order, Discrepancy) %>%
  mutate(row_id = paste(Order, Discrepancy, sep = "__"))

all_df$row_id <- factor(paste(all_df$Order, all_df$Discrepancy, sep = "__"),
                         levels = rev(row_order$row_id))

# Matches the manuscript's own legend badges: NPP = black (the
# baseline vline below), NPP-SEQ = pink/magenta, ONPP = blue,
# ONPP-SEQ = orange.
method_colors <- c("NPP_SEQ" = "#C04090", "ONPP" = "#0072B2", "ONPP_SEQ" = "#D95F02")
method_labels <- c("NPP_SEQ" = "NPP-SEQ", "ONPP" = "ONPP", "ONPP_SEQ" = "ONPP-SEQ")

#-------------------------------------------------------------------------------
# Alternating background bands (one per Order group, spanning its 2-3
# sub-rows) plus a highlighted band on "Incompatible, none" -- the
# scenario class that breaks ONPP's order restriction, singled out the
# same way elsewhere in this project (prior_comparison_*.R's header
# comments).
#-------------------------------------------------------------------------------

order_band_df <- row_order %>%
  mutate(y_num = as.numeric(factor(row_id, levels = rev(row_order$row_id)))) %>%
  group_by(Order) %>%
  summarise(ymin = min(y_num) - 0.5, ymax = max(y_num) + 0.5, .groups = "drop") %>%
  mutate(fill_group = as.numeric(factor(Order, levels = order_levels)) %% 2)

special_band_df <- row_order %>%
  filter(Order == "Incompatible", Discrepancy == "none") %>%
  mutate(y_num = as.numeric(factor(row_id, levels = rev(row_order$row_id))),
         ymin = y_num - 0.5, ymax = y_num + 0.5)

# Thin separator line between every consecutive pair of scenario rows
# (not just between Order groups, which the bands above already set off).
scenario_sep_y <- seq(1.5, nrow(row_order) - 0.5, by = 1)

# Bold separator line between consecutive ORDER groups (Compatible /
# Semi-compat. / Incompatible / Neutral / Cancellation) -- order_band_df's
# rows are already in order_levels order (group_by/summarise on a factor
# sorts by level), so each group's own ymax is the boundary below it;
# drop the last one, which is the outer edge rather than a boundary.
order_sep_y <- order_band_df$ymax[-1]

# Built directly from all_df's own NA ratios (not a separate cross-join),
# so this is guaranteed to land only on the panels/rows that are actually
# NA -- the Regression panels' Cancellation rows -- and never bleed into
# Bernoulli/Normal, which have real Cancellation data.
not_defined_df <- all_df %>%
  filter(is.na(ratio)) %>%
  distinct(row_id, panel)

#-------------------------------------------------------------------------------
# Main 4-panel plot.
#-------------------------------------------------------------------------------

p_main <- ggplot(all_df, aes(x = ratio, y = row_id)) +
  geom_rect(data = order_band_df,
            aes(xmin = -Inf, xmax = Inf, ymin = ymin, ymax = ymax, fill = factor(fill_group)),
            inherit.aes = FALSE, alpha = 0.6) +
  geom_rect(data = special_band_df,
            aes(xmin = -Inf, xmax = Inf, ymin = ymin, ymax = ymax),
            inherit.aes = FALSE, fill = "#F7C6CB", alpha = 0.7) +
  scale_fill_manual(values = c("0" = "white", "1" = "grey93"), guide = "none") +
  geom_hline(yintercept = scenario_sep_y, color = "grey80", linewidth = 0.3) +
  geom_hline(yintercept = order_sep_y, color = "black", linewidth = 0.9) +
  geom_vline(xintercept = 1, color = "black", linewidth = 0.7) +
  geom_point(aes(color = method), size = 3.4, na.rm = TRUE) +
  geom_text(data = not_defined_df, aes(x = 1, y = row_id, label = "not defined"),
            inherit.aes = FALSE, color = "grey55", fontface = "italic", size = 4) +
  scale_color_manual(values = method_colors, labels = method_labels, name = NULL) +
  scale_x_log10(breaks = c(0.5, 0.67, 1, 1.5, 2), labels = c("0.5", "0.67", "1", "1.5", "2")) +
  scale_y_discrete(labels = function(x) sub("^.*__", "", x)) +
  facet_wrap(~ panel, nrow = 1) +
  labs(x = "MSE / MSE of NPP", y = NULL) +
  .sce3_theme() +
  theme(
    legend.position    = "none",
    strip.text         = element_text(face = "bold", size = 15, hjust = 0),
    panel.grid.major.y = element_blank(),
    axis.text.y        = element_text(size = 12),
    axis.text.x        = element_text(size = 11),
    axis.title          = element_text(size = 14)
  )

#-------------------------------------------------------------------------------
# Narrow left-hand strip carrying the Order group names ONCE (not repeated
# per panel) -- plain text at each group's vertical midpoint, no axis/
# border, glued to p_main with patchwork so its rows line up exactly.
#-------------------------------------------------------------------------------

order_label_df <- order_band_df %>% mutate(y_mid = (ymin + ymax) / 2)

p_order_labels <- ggplot(order_label_df, aes(x = 0, y = y_mid, label = Order)) +
  geom_text(hjust = 1, size = 4.6, color = "#2C5F8A", fontface = "bold") +
  coord_cartesian(xlim = c(-1, 0), ylim = c(0.5, nrow(row_order) + 0.5), clip = "off") +
  theme_void()

p_combined <- p_order_labels + p_main + plot_layout(widths = c(0.14, 1))

print(p_combined)
ggsave("results/figures/poster/combined_mse_ratio_forest_plot.pdf",
       p_combined, width = 13, height = 4, device = cairo_pdf)

# ---- To add specific result callouts once you've inspected the actual
# numbers (e.g. "SEQ prior pulls sigma^2 toward its prior mean -- backfires
# under Cancellation", as in the reference figure), add annotate() calls
# with real coordinates read off `all_df` after running this script --
# left as a placeholder here since those are findings about THIS run's
# results, not something to hardcode sight unseen. ----
