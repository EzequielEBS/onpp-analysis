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
library(scales)
library(RColorBrewer)

# run auxiliary functions
source("code/aux_fun_sim.R")

gamma_model_lm <- cmdstan_model("code/models/lm/gamma_lm.stan")
delta_model_lm <- cmdstan_model("code/models/lm/eta_lm.stan")

load("../../simulated-data-regression/data/true_params.RData")
betastar <- beta[[1]]
sgstar <- 1

r <- 3

# The 12 sce1.1..sce4.3 objects below used to be ~170 lines of copy-pasted
# `list(par = list(p=3, a=2, ..., post=0), data = qs_read(...)$data[r:(r+1)])`
# blocks -- identical `par` in every block, differing only in which
# sce<Order>_<Congruence>_varying_beta_n0_ge_n.qs2 file gets read. Built
# through the same order/congruence loop used in prepare_data_lm.r /
# simulations_lm.r / analysis_lm.R instead. Also fixes a small inefficiency:
# `$data[r:(r+1)]` sliced out *two* replicates (r and r+1) as a length-2
# list even though only the first ([[1]]) was ever used downstream (see
# sample_prior_post_lm() below, which now takes a single replicate
# directly) -- `$data[[r]]` pulls just the one replicate actually needed.
lm_prior_par <- list(p = 3, a = 2, b = 1, V0 = 100*diag(3), mu0 = rep(0, 3),
                      al = 1, bl = 1, alpha = rep(1, 4), post = 0)

sce_list <- list()
for (order_i in 1:4) {
  for (cong_j in 1:3) {
    order_roman <- as.character(as.roman(order_i))
    cong_roman  <- as.character(as.roman(cong_j))
    sce_file <- sprintf("results/sim_data/lm/sce%s_%s_varying_beta_n0_ge_n.qs2", order_roman, cong_roman)
    sce_list[[length(sce_list) + 1]] <- list(
      par = lm_prior_par,
      data = qs_read(sce_file)$data[[r]]
    )
  }
}

# Loads the 12 sim<order>_<cong>_lm_varying_beta_n0_ge_n.qs2 files (order
# I-IV outer, congruence I-III inner -- same convention as analysis_lm.R's
# load_lm_sims()), replacing another 12 copy-pasted qs_read() calls.
sim_sces_n0_ge_n <- list()
for (order_i in 1:4) {
  for (cong_j in 1:3) {
    file <- sprintf("results/samples/lm/sim%d_%d_lm_varying_beta_n0_ge_n.qs2", order_i, cong_j)
    sim_sces_n0_ge_n[[length(sim_sces_n0_ge_n) + 1]] <- qs_read(file)
  }
}

df_sce <- data.frame(
  Eps = c(
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
  ),
  Order       = c(rep("Compatible",    3),
                      rep("Semi-compat.",  3),
                      rep("Incompatible",  3),
                      rep("Neutral",       3)
                    ),
  Discrepancy = rep(c("None","Small","Large"), 4)
                  )


sample_prior_post_lm <- function(i, r, sce_list, sim_list) {
  sce <- sce_list[[i]]
  sim <- sim_list[[i]]
  # This only ever needs ONE prior draw, but used to get it via
  # sim_sce(model="lm", 1, 1, ...) -- which spins up a 1-worker PSOCK
  # cluster and reloads cmdstanr/dplyr/posterior on it just to run a single
  # Stan fit. sample_sce_lm() is the same per-replicate fit sim_sce() calls
  # internally; calling it directly (sce$par$post is already forced to 0
  # above, so this draws from the prior) gives the same single draw without
  # any cluster setup/teardown overhead.
  #
  # Cached to disk (qs2) so re-running this script -- e.g. just to tweak
  # figure aesthetics -- doesn't re-fit the Stan prior model every time.
  # Keyed on r too (unlike bin/normal) since lm's sce_list$data is itself
  # built from replicate r (see sce_list construction above). Delete the
  # cache file to force a re-sample (e.g. after changing a scenario's
  # parameters or r).
  cache_file <- sprintf("results/samples/lm/prior_cache_sce%d_r%d.qs2", i, r)
  if (file.exists(cache_file)) {
    prior_fit <- qs2::qs_read(cache_file)
  } else {
    prior_fit <- sample_sce_lm(sce$par, sce$data, gamma_model_lm, delta_model_lm)
    qs2::qs_save(prior_fit, file = cache_file)
  }
  eta_prior <- prior_fit[c("delta_npp", "delta_nppseq", "delta_onpp", "delta_onppseq")]
  theta_prior <- prior_fit[c("theta_npp", "theta_nppseq", "theta_onpp", "theta_onppseq")]
  eta_post <- sim$delta[[r]]
  theta_post <- sim$theta[[r]]
  list(
    eta_prior = eta_prior,
    theta_prior = theta_prior,
    eta_post = eta_post,
    theta_post = theta_post
  )
}

# OLS/MLE point estimates of (beta0, beta1, beta2) and the MLE residual
# variance, computed separately from each historical arm and from the
# current data, for marking on the theta prior/posterior plots below.
# `data` is one replicate's entry from sce_list[[i]]$data (as produced by
# prepare_data_stan() in prepare_data_lm.r): the historical arms are stored
# flattened (X0_flat/y0_flat + startid_*/dims_*) since that's the shape
# Stan needs, so they're unflattened back into per-arm design matrices
# here; the current arm is already a plain n x p matrix / n x 1 matrix.
compute_lm_mle <- function(data) {
  unflatten <- function(flat, startid, dims, k) {
    len <- dims[k, 1] * dims[k, 2]
    matrix(flat[startid[k]:(startid[k] + len - 1)], nrow = dims[k, 1], ncol = dims[k, 2])
  }
  ols <- function(X, y) {
    fit <- lm.fit(X, y)
    list(beta = coef(fit), sigma2 = mean(resid(fit)^2))
  }

  K <- nrow(data$dims_X0)
  hist_est <- lapply(seq_len(K), function(k) {
    X0_k <- unflatten(data$X0_flat, data$startid_X0, data$dims_X0, k)
    y0_k <- as.vector(unflatten(data$y0_flat, data$startid_y0, data$dims_y0, k))
    ols(X0_k, y0_k)
  })
  cur_est <- ols(as.matrix(data$X), as.vector(data$y))

  est <- c(hist_est, list(cur_est))
  data.frame(
    source = c(paste0("Hist ", seq_len(K)), "Current"),
    beta0  = vapply(est, function(e) e$beta[1], numeric(1)),
    beta1  = vapply(est, function(e) e$beta[2], numeric(1)),
    beta2  = vapply(est, function(e) e$beta[3], numeric(1)),
    sigma2 = vapply(est, function(e) e$sigma2, numeric(1))
  )
}

plot_theta_lm <- function(
  draws,
  ordered = F,
  post = F,
  theta_true = c(betastar, sgstar),
  upp_band,
  low_band,
  mle = NULL,
  title = NULL,
  leg = T
) {
  # Colors for the SEQ/Default density curves, extended with one more color
  # per historical arm plus current data when `mle` marks are requested.
  # SEQ/Default now use the same black/reddish-purple pair as
  # .default_seq_colors in aux_fun_sim.R (plot_theta_prior_post()/
  # plot_eta_prior_post(), bin and normal's shared scalar-theta/eta plots),
  # instead of the plain "blue"/"red" names this file used before -- keeps
  # the "Default vs SEQ" color coding consistent across all three families'
  # figures, and deliberately distinct from plot_sce3_theta()'s NPP/ONPP
  # blue/vermillion pair (same color meaning two different comparisons
  # across figures was confusing). The 4 MLE-marker colors are unrelated to
  # that pair and were already chosen (dataviz skill categorical slots) and
  # validated colorblind-safe as a 6-color set with
  # scripts/validate_palette.js --pairs all (worst normal-vision pair
  # #008300<->#1baf7a at ~15.6, just above the >=15 floor; the one pair in
  # the 6-8 CVD band, aqua/magenta, is covered by the direct legend labels
  # below, which is the documented mitigation for that band). Assumes at
  # most 3 historical arms (K = 3 throughout this file's scenarios).
  #
  # NOTE: the inline color comments below previously said "purple" for
  # #00FF00 and "orange" for #00FFFF -- neither matches those hex values
  # (pure green and cyan, respectively). Corrected; the actual colors used
  # were always correct, only the comments describing them were wrong.
  color_values <- c("SEQ" = unname(.default_seq_colors["SEQ"]), "pi" = unname(.default_seq_colors["Default"]))
  color_labels <- c("SEQ" = "SEQ", "pi" = "Default")

  if (!is.null(mle)) {
    color_values <- c(color_values,
                      "Hist 1" = "#FFFF00",
                      "Hist 2" = "#00FF00",  # green
                      "Hist 3" = "#00FFFF",  # cyan
                      "Current" = "#FF00FF") # pink/magenta
    color_labels <- c(color_labels,
                      "Hist 1" = "Hist. 1 MLE", "Hist 2" = "Hist. 2 MLE",
                      "Hist 3" = "Hist. 3 MLE", "Current" = "Current MLE")
  }
  draws_beta_nppseq <- draws$theta_nppseq[, 1:3]
  draws_beta_npp <- draws$theta_npp[, 1:3]
  draws_s2_nppseq <- draws$theta_nppseq[, 4]
  draws_s2_npp <- draws$theta_npp[, 4]

  if (ordered) {
    draws_beta_nppseq <- draws$theta_onppseq[, 1:3]
    draws_beta_npp <- draws$theta_onpp[, 1:3]
    draws_s2_nppseq <- draws$theta_onppseq[, 4]
    draws_s2_npp <- draws$theta_onpp[, 4]
  }

  if (post) {
    y_lab <- expression(paste(pi[0], "(", theta ~ "|" ~ D ~ "," ~ D[0], ")"))
  } else {
    y_lab <- expression(paste(pi[0], "(", theta, ")"))
  }

  plot_beta1 <- ggplot() +
    # Prior curves dotted, posterior curves solid -- same reasoning as
    # aux_fun_sim.R's plot_eta_prior_post()/plot_theta_prior_post().
    geom_density(aes(x = draws_beta_nppseq[, 1], color = "SEQ"), linewidth = 0.8,
                linetype = if (post) "solid" else "dashed") +
    geom_density(aes(x = draws_beta_npp[, 1], color = "pi"), linewidth = 0.8,
                linetype = if (post) "solid" else "dashed") +
    labs(x = expression(beta[0]), y = y_lab) +
    # grey, not black, so it doesn't visually collide with the (now black) "Default" curve
    geom_vline(xintercept = theta_true[1], linetype = "dotted", color = "grey40") +
    scale_color_manual(name = NULL, values = color_values,
                        labels = color_labels) +
    scale_y_continuous(labels = function(x) format(round(x, 1), nsmall = 1)) +
    scale_x_continuous(labels = function(x) format(round(x, 1), nsmall = 1),
                        limits = c(theta_true[1] - low_band[1], theta_true[1] + upp_band[1])
                      ) +
    .sce3_theme() +
    ggtitle(
      title
    )
  if (!is.null(mle)) {
    plot_beta1 <- plot_beta1 +
      geom_point(data = mle, aes(x = beta0, y = 0, color = source),
                 shape = 18, size = 4, stroke = 1.2)
  }
  plot_beta2 <- ggplot() +
    # Prior curves dotted, posterior curves solid -- same reasoning as
    # aux_fun_sim.R's plot_eta_prior_post()/plot_theta_prior_post().
    geom_density(aes(x = draws_beta_nppseq[, 2], color = "SEQ"), linewidth = 0.8,
                linetype = if (post) "solid" else "dashed") +
    geom_density(aes(x = draws_beta_npp[, 2], color = "pi"), linewidth = 0.8,
                linetype = if (post) "solid" else "dashed") +
    labs(x = expression(beta[1]), y = "") +
    geom_vline(xintercept = theta_true[2], linetype = "dotted", color = "grey40") +
    scale_color_manual(name = NULL, values = color_values,
                        labels = color_labels) +
    scale_y_continuous(labels = function(x) format(round(x, 1), nsmall = 1)) +
    scale_x_continuous(labels = function(x) format(round(x, 1), nsmall = 1),
                        limits = c(theta_true[2] - low_band[2], theta_true[2] + upp_band[2])
                      ) +
    .sce3_theme()
  if (!is.null(mle)) {
    plot_beta2 <- plot_beta2 +
      geom_point(data = mle, aes(x = beta1, y = 0, color = source),
                 shape = 18, size = 4, stroke = 1.2)
  }
  plot_beta3 <- ggplot() +
    # Prior curves dotted, posterior curves solid -- same reasoning as
    # aux_fun_sim.R's plot_eta_prior_post()/plot_theta_prior_post().
    geom_density(aes(x = draws_beta_nppseq[, 3], color = "SEQ"), linewidth = 0.8,
                linetype = if (post) "solid" else "dashed") +
    geom_density(aes(x = draws_beta_npp[, 3], color = "pi"), linewidth = 0.8,
                linetype = if (post) "solid" else "dashed") +
    geom_vline(xintercept = theta_true[3], linetype = "dotted", color = "grey40") +
    labs(x = expression(beta[2]), y = "") +
    scale_color_manual(name = NULL, values = color_values,
                        labels = color_labels) +
    scale_y_continuous(labels = function(x) format(round(x, 1), nsmall = 1)) +
    scale_x_continuous(labels = function(x) format(round(x, 1), nsmall = 1),
                        limits = c(theta_true[3] - low_band[3], theta_true[3] + upp_band[3])
                      ) +
    .sce3_theme()
  if (!is.null(mle)) {
    plot_beta3 <- plot_beta3 +
      geom_point(data = mle, aes(x = beta2, y = 0, color = source),
                 shape = 18, size = 4, stroke = 1.2)
  }
  plot_s2 <- ggplot() +
    # Prior curves dotted, posterior curves solid -- same reasoning as
    # aux_fun_sim.R's plot_eta_prior_post()/plot_theta_prior_post().
    geom_density(aes(x = draws_s2_nppseq, color = "SEQ"), linewidth = 0.8,
                linetype = if (post) "solid" else "dashed") +
    geom_density(aes(x = draws_s2_npp, color = "pi"), linewidth = 0.8,
                linetype = if (post) "solid" else "dashed") +
    geom_vline(xintercept = theta_true[4], linetype = "dotted", color = "grey40") +
    labs(x = expression(sigma^2), y = "") +
    scale_color_manual(name = NULL, values = color_values,
                        labels = color_labels) +
    scale_y_continuous(labels = function(x) format(round(x, 1), nsmall = 1)) +
    scale_x_continuous(labels = function(x) format(round(x, 1), nsmall = 1),
                        limits = c(theta_true[4] - low_band[4], theta_true[4] + upp_band[4])
                      ) +
    .sce3_theme()
  if (!is.null(mle)) {
    plot_s2 <- plot_s2 +
      geom_point(data = mle, aes(x = sigma2, y = 0, color = source),
                 shape = 18, size = 4, stroke = 1.2)
  }

  out <- plot_beta1 + plot_beta2 + plot_beta3 + plot_s2 +
  plot_layout(ncol = 4, guides = "collect")
  if (leg) {
    out <- out & theme(legend.position = "top")
  } else {
    out <- out & theme(legend.position = "none")
  }
  out
}

# plot_eta_lm() used to be defined here (byte-for-byte identical, apart
# from the family name, to analysis_seq_bin.r's and
# analysis_seq_normal.r's own copies). It's now plot_eta_prior_post() in
# aux_fun_sim.R -- see the comment above that function for the full
# rationale. plot_theta_lm() below stays local: lm's theta is a 4-vector
# (beta0/beta1/beta2/sigma^2) with per-parameter asymmetric bands and an
# optional MLE-marker overlay, different enough from bin/normal's scalar
# theta that it isn't a good fit for the shared plot_theta_prior_post().

plots_npp_sce <- lapply(1:3, function(i) {
  samples <- sample_prior_post_lm(i, r, sce_list, sim_sces_n0_ge_n)
  # MLE marks (per historical arm + current data) for the theta plots below
  # -- computed from the same replicate's data that fits sce_list[[i]] used,
  # so the marks stay correct even for the reordered-arm (Order II/III)
  # scenarios.
  mle_i <- compute_lm_mle(sce_list[[i]]$data)
  # No "Scenario ..." label -- epsilon/Order/Discrepancy already identify
  # the scenario.
  title <- bquote(
    epsilon == .(df_sce$Eps[i])~", Order:"~.(df_sce$Order[i])~", Discrepancy:"~.(df_sce$Discrepancy[i])
  )
  if (i == 1){
    plot_theta_prior <- plot_theta_lm(samples$theta_prior, title = title, leg = T, mle = mle_i,
      upp_band = c(15,19,25,2), low_band = c(24,24,22,1)
    )
    plot_theta_post <- plot_theta_lm(samples$theta_post, post = T, leg = F, mle = mle_i,
      upp_band = c(0.75,0.75,0.75,2), low_band = c(0.75,0.75,0.75,1)
    )
  } else {
    plot_theta_prior <- plot_theta_lm(samples$theta_prior, title = title, leg = F, mle = mle_i,
      upp_band = c(15,19,25,2), low_band = c(24,24,22,1)
    )
    plot_theta_post <- plot_theta_lm(samples$theta_post, post = T, leg = F, mle = mle_i,
      upp_band = c(0.75,0.75,0.75,2), low_band = c(0.75,0.75,0.75,1)
    )
  }

  if (i != 3){
    plot_eta_prior <- plot_eta_prior_post(samples$eta_prior, title = title, leg = F)
    plot_eta_post <- plot_eta_prior_post(samples$eta_post, post = T, leg = F)
  } else {
    plot_eta_prior <- plot_eta_prior_post(samples$eta_prior, title = title, leg = F)
    plot_eta_post <- plot_eta_prior_post(samples$eta_post, post = T, leg = T)
  }
  
  plot_theta <- plot_theta_prior / plot_theta_post
  # No plot_layout(guides = "collect") here -- see the note in
  # analysis_seq_bin.r on why relying on collect produced duplicate/
  # misplaced legends.
  plot_eta <- plot_eta_prior / plot_eta_post
  list(
    plot_theta = plot_theta,
    plot_eta = plot_eta
    )
})

# No guides = "collect" (or `&` position override) here -- see the full
# note in analysis_seq_bin.r. Only scenario 3's post plot sets leg = T,
# so it's the only one that builds a legend, and it already renders at
# the bottom of the whole figure since it's the last block.
plots_eta_npp <- lapply(1:3, function(i) plots_npp_sce[[i]]$plot_eta) %>%
  wrap_plots(ncol = 1)
plots_theta_npp <- lapply(1:3, function(i) plots_npp_sce[[i]]$plot_theta) %>%
  wrap_plots(ncol = 1)

ggsave("results/figures/lm/theta_prior_post_npp.pdf", plots_theta_npp, width = 12, height = 13.5)
ggsave("results/figures/lm/eta_prior_post_npp.pdf", plots_eta_npp, width = 9, height = 13.5)

#------------------------------------------------------------------------------
# Small sample size
#------------------------------------------------------------------------------
# NOTE: this section header was already here with nothing under it -- no
# code follows in the original file either. Left as-is (not implemented
# here) rather than guessing at what the small-size prior/posterior
# comparison was meant to look like; the small-size scenario data/samples
# already exist (sceI_I/III_..._small_size.qs2, sim1_1/3_..._small_size.qs2
# from prepare_data_lm.r/simulations_lm.r) if this gets picked up later.
