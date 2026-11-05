# ==============================================================================
# Posterior density plots (eta + theta) for one replicate, cancellation-
# effect scenarios, combined across models
# ==============================================================================
#
# The "cancellation effect" scenarios (theta0 settings V.I/V.II -- the ones
# that don't fit the order x congruence grid every other scenario uses) only
# exist for the bernoulli (bin) and normal-fixed-variance (normal) families;
# lm has no equivalent (its scenario table only covers the 4 order x 3
# congruence settings -- see simulations_lm.r). So there's nothing to plot
# for lm here.
#
# Reuses plot_eta_theta_density_styled() (defined in aux_fun_sim.R), which
# builds a single combined figure for one simulation replicate -- an eta
# density panel faceted by method (NPP/NPP-SEQ/ONPP/ONPP-SEQ) on top, theta
# density below -- colored by eta component / method respectively, in the
# same visual language (.sce3_theme(), Okabe-Ito colorblind palette) as this
# file's other sce3/prior-post figures.
#
# For each of the two sub-scenarios (V.I, V.II), the bin and normal figures
# are stacked into one combined PDF (bin on top, normal below) via
# patchwork's `/`, rather than saved as four separate files -- same theta0
# pattern (one setting pushed to an extreme, one pulled back to cancel it
# out), so seeing both models' behavior on that pattern side by side in one
# file is more useful than flipping between two.

library(cmdstanr)
library(dplyr)
library(tidyverse)
library(posterior)
library(patchwork)
library(qs2)

source("code/aux_fun_sim.R")

replicate <- 1  # which simulated replicate to plot

# ------------------------------------------------------------------------------
# Bernoulli (bin) -- theta0 = (0.10, 0.90, 0.50) / (0.50, 0.90, 0.10),
# true theta = 0.5, both "Large" discrepancy (scenarios_bin.R / df_mse_theta
# rows 13-14 in analysis_bin.r)
# ------------------------------------------------------------------------------
sim5_1_bin <- qs2::qs_read("results/samples/bin/sim5_1_bin.qs2")
sim5_2_bin <- qs2::qs_read("results/samples/bin/sim5_2_bin.qs2")

title_canc_I_bin <- bquote(
  "Bernoulli, " * theta[0] == "(0.10, 0.90, 0.50)" * ", Order: Cancellation effect, Discrepancy: Large"
)
plot_canc_I_bin <- plot_eta_theta_density_styled(
  sim5_1_bin, replicate = replicate, family = "bin", true_theta = 0.5,
  title = title_canc_I_bin
)

title_canc_II_bin <- bquote(
  "Bernoulli, " * theta[0] == "(0.50, 0.90, 0.10)" * ", Order: Cancellation effect, Discrepancy: Large"
)
plot_canc_II_bin <- plot_eta_theta_density_styled(
  sim5_2_bin, replicate = replicate, family = "bin", true_theta = 0.5,
  title = title_canc_II_bin
)

# ------------------------------------------------------------------------------
# Normal (fixed variance) -- theta0 = (3, 7, 5) / (5, 7, 3), true theta = 5
# (simulations_normal.r's normal_default_par), both "Large" discrepancy
# ------------------------------------------------------------------------------
sim5_1_normal <- qs2::qs_read("results/samples/normal/sim5_1_normal_fixed_var.qs2")
sim5_2_normal <- qs2::qs_read("results/samples/normal/sim5_2_normal_fixed_var.qs2")

title_canc_I_normal <- bquote(
  "Normal, " * theta[0] == "(3.0, 7.0, 5.0)" * ", Order: Cancellation effect, Discrepancy: Large"
)
plot_canc_I_normal <- plot_eta_theta_density_styled(
  sim5_1_normal, replicate = replicate, family = "normal", true_theta = 5,
  title = title_canc_I_normal
)

title_canc_II_normal <- bquote(
  "Normal, " * theta[0] == "(5.0, 7.0, 3.0)" * ", Order: Cancellation effect, Discrepancy: Large"
)
plot_canc_II_normal <- plot_eta_theta_density_styled(
  sim5_2_normal, replicate = replicate, family = "normal", true_theta = 5,
  title = title_canc_II_normal
)

# ------------------------------------------------------------------------------
# Combine bin + normal into one figure per sub-scenario
# ------------------------------------------------------------------------------
# Each plot_eta_theta_density_styled() call already bakes in its own
# plot_layout(guides = "collect") legend, so stacking two of them with a
# bare / lets the two blocks end up different sizes if their internal
# legend space isn't identical (this is what made the bin block visibly
# taller than the normal block). Using "/" + plot_layout(heights = ...)
# re-flows that heights spec into the *nested* layouts inside each
# already-assembled block too (each of which has its own multi-panel
# design internally), which is what caused wrap_dims() to error looking
# for a 3rd panel that doesn't exist. wrap_plots() with ncol = 1 instead
# treats plot_canc_*_bin/normal as two single opaque units and only
# applies heights = c(1, 1) to that outer 2-row layout, leaving each
# block's own internal layout untouched.
plot_canc_I  <- wrap_plots(plot_canc_I_bin,  plot_canc_I_normal,  ncol = 1, heights = c(1, 1))
plot_canc_II <- wrap_plots(plot_canc_II_bin, plot_canc_II_normal, ncol = 1, heights = c(1, 1))

ggsave("results/figures/cancellation_density_sceV_I_combined.pdf",
       plot_canc_I, width = 10, height = 16)
ggsave("results/figures/cancellation_density_sceV_II_combined.pdf",
       plot_canc_II, width = 10, height = 16)
