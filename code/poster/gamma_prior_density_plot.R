# =============================================================================
# Marginal prior densities of eta_1, eta_2, eta_3 induced by two choices of
# Dirichlet concentration on gamma -- the "couple AND order" ONPP prior
# (gamma_bin.stan / gamma_lm.stan: gamma ~ Dirichlet(alpha), eta =
# cumulative_sum(gamma)[1:3], with eta_4 = 1 implicit and dropped) -- for
# K = 3 historical datasets, so gamma is 4-dimensional.
#
# Two concentration vectors compared, both symmetric across the 4
# components:
#   - "Han et al.": alpha = rep(1/4, 4) -- a symmetric Dirichlet with
#     concentration < 1, which pushes mass toward the simplex's corners/
#     edges (U-shaped marginals, blowing up at 0/1).
#   - "flat": alpha = rep(1, 4) -- the project's own ONPP default (see
#     prior_comparison_sample.R's `alpha = rep(1, 4)`), the uniform
#     distribution on the simplex (unimodal, bell-shaped marginals).
#
# No simulation needed: for gamma ~ Dirichlet(alpha_1, ..., alpha_K), the
# partial sum eta_k = gamma_1 + ... + gamma_k is EXACTLY
#   eta_k ~ Beta(sum(alpha[1:k]), sum(alpha[(k+1):K]))
# (a standard Dirichlet/Beta aggregation identity), so each panel below is
# just three dbeta() curves evaluated on a fine grid -- this is the prior
# ITSELF, not a Monte Carlo approximation of it.
# =============================================================================

library(dplyr)
library(tidyverse)
library(patchwork)

# Reuses .sce3_theme() for the shared visual language (matches every other
# sce3/prior figure in this project).
source("code/aux_fun_sim.R")

#-------------------------------------------------------------------------------
# Marginal eta_1..eta_{K-1} densities for gamma ~ Dirichlet(alpha), via the
# partial-sum-of-Dirichlet-is-Beta identity above. Returns one row per
# (eta grid point, component).
#-------------------------------------------------------------------------------

eta_marginal_density_df <- function(alpha, n_grid = 1000) {
  K <- length(alpha)
  total <- sum(alpha)
  eta_grid <- seq(0.0001, 0.9999, length.out = n_grid)

  do.call(rbind, lapply(seq_len(K - 1), function(k) {
    shape1 <- sum(alpha[seq_len(k)])
    shape2 <- total - shape1
    data.frame(
      eta       = eta_grid,
      density   = dbeta(eta_grid, shape1, shape2),
      component = k
    )
  }))
}

K <- 3  # number of historical datasets -- gamma is (K+1)-dimensional

prior_specs <- list(
  "Han et al." = list(alpha = rep(1 / 4, K + 1),
                       title = bquote(gamma ~ "~" ~ Dir(textstyle(frac(1, 4)) * "," ~ "..." * "," ~ textstyle(frac(1, 4))) ~ "(Han et al.)")),
  "Flat"       = list(alpha = rep(1, K + 1),
                       title = bquote(gamma ~ "~" ~ Dir(1 * "," ~ 1 * "," ~ 1 * "," ~ 1) ~ "(flat)"))
)

eta_long <- dplyr::bind_rows(lapply(names(prior_specs), function(nm) {
  df <- eta_marginal_density_df(prior_specs[[nm]]$alpha)
  df$prior <- nm
  df
}))
eta_long$prior <- factor(eta_long$prior, levels = names(prior_specs))
eta_long$component <- factor(eta_long$component, levels = seq_len(K))

#-------------------------------------------------------------------------------
# Blue gradient (light -> dark navy, eta_1 -> eta_3), matching the project's
# convention of a sequential palette for an ordered set of components.
#-------------------------------------------------------------------------------

eta_colors <- c("1" =  "#E69F00", "2" = "#56B4E9", "3" = "#009E73")
eta_labels <- c("1" = "eta[1]", "2" = "eta[2]", "3" = "eta[3]")

# Direct curve labels (component name placed right on/above its own curve)
# instead of a legend -- matches the reference figure. Placement is
# hand-picked per prior x component since the two priors' curves have very
# different shapes (U-shaped vs bell-shaped).
label_positions <- dplyr::bind_rows(
  data.frame(prior = "Han et al.", component = factor(c(1, 2, 3)),
             eta = c(0.14, 0.50, 0.86),
             density = c(dbeta(0.14, 0.25, 0.75) + 0.3,
                         dbeta(0.50, 0.50, 0.50) + 0.3,
                         dbeta(0.86, 0.75, 0.25) + 0.3)),
  data.frame(prior = "Flat", component = factor(c(1, 2, 3)),
             eta = c(0.08, 0.50, 0.92),
             density = c(dbeta(0.08, 1, 3) + 0.15,
                         dbeta(0.50, 2, 2) + 0.15,
                         dbeta(0.92, 3, 1) + 0.15))
)
label_positions$prior <- factor(label_positions$prior, levels = names(prior_specs))

# `ylim_cap`: Beta(alpha, beta) with alpha or beta < 1 (the "Han et al."
# case, shapes 0.25/0.75) diverges to infinity at the boundary it's < 1
# on -- the auto y-range needed to fit that spike (order ~200) crushes
# every curve flat near y = 0 everywhere else, which is what made the
# three curves look "together"/indistinguishable in the middle. Capping
# the y-axis with coord_cartesian() (not scale_y_continuous(limits=...),
# which would DROP the out-of-range points and break the lines) lets the
# spikes run off the top of the panel -- same look as the reference
# figure -- while the mid-range, where the three curves actually differ
# in shape, is no longer squashed into a sliver at the bottom.
make_panel <- function(prior_name, ylim_cap = NA) {
  df_panel    <- eta_long[eta_long$prior == prior_name, ]
  labels_panel <- label_positions[label_positions$prior == prior_name, ]

  ggplot(df_panel, aes(x = eta, y = density, color = component)) +
    geom_line(linewidth = 1) +
    geom_text(data = labels_panel,
              aes(x = eta, y = density, color = component,
                  label = eta_labels[as.character(component)]),
              parse = TRUE, size = 4, fontface = "italic", show.legend = FALSE) +
    scale_color_manual(values = eta_colors, guide = "none") +
    coord_cartesian(ylim = c(0, ylim_cap)) +
    labs(x = expression(eta[k]), y = NULL,
         title = prior_specs[[prior_name]]$title) +
    .sce3_theme() +
    theme(plot.title = element_text(size = 11),
          axis.text  = element_text(size = 12),
          axis.title = element_text(size = 13))
}

p_gamma_prior <- make_panel("Han et al.", ylim_cap = 6) + make_panel("Flat")

print(p_gamma_prior)
ggsave("results/figures/poster/gamma_prior_eta_density.pdf", p_gamma_prior, width = 10, height = 3)
