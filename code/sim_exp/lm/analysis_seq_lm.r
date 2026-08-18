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
lm_prior_par <- list(p = 3, a = 2, b = 1, V0 = diag(3), mu0 = rep(0, 3),
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
  prior_fit <- sample_sce_lm(sce$par, sce$data, gamma_model_lm, delta_model_lm)
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

plot_theta_lm <- function(
  draws, 
  ordered = F,
  post = F,
  theta_true = c(betastar, sgstar), 
  # c = c(3, 3), 
  title = NULL,
  leg = T
) {
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
    geom_density(aes(x = draws_beta_nppseq[, 1], color = "SEQ"), size = 1) +
    geom_density(aes(x = draws_beta_npp[, 1], color = "pi"), size = 1) +
    labs(x = expression(beta[0]), y = y_lab) +
    geom_vline(xintercept = theta_true[1], linetype = "dashed", color = "black") +
    scale_color_manual(name = NULL, values = c("SEQ" = "blue", "pi" = "red"), 
                        labels = c("SEQ" = "SEQ",
                                    "pi" = "Default")) +
    scale_y_continuous(labels = function(x) format(round(x, 1), nsmall = 1)) +
    scale_x_continuous(labels = function(x) format(round(x, 1), nsmall = 1)
                        # limits = c(theta_true[1] - c[1], theta_true[1] + c[1])
                      ) +
    theme_minimal() +
    ggtitle(
      title
    )
  plot_beta2 <- ggplot() +
    geom_density(aes(x = draws_beta_nppseq[, 2], color = "SEQ"), size = 1) +
    geom_density(aes(x = draws_beta_npp[, 2], color = "pi"), size = 1) +
    labs(x = expression(beta[1]), y = "") +
    geom_vline(xintercept = theta_true[2], linetype = "dashed", color = "black") +
    scale_color_manual(name = NULL, values = c("SEQ" = "blue", "pi" = "red"), 
                        labels = c("SEQ" = "SEQ",
                                    "pi" = "Default")) +
    scale_y_continuous(labels = function(x) format(round(x, 1), nsmall = 1)) +
    scale_x_continuous(labels = function(x) format(round(x, 1), nsmall = 1)
                        # limits = c(theta_true[2] - c[1], theta_true[2] + c[1])
                      ) +
    theme_minimal()
  plot_beta3 <- ggplot() +
    geom_density(aes(x = draws_beta_nppseq[, 3], color = "SEQ"), size = 1) +
    geom_density(aes(x = draws_beta_npp[, 3], color = "pi"), size = 1) +
    geom_vline(xintercept = theta_true[3], linetype = "dashed", color = "black") +
    labs(x = expression(beta[2]), y = "") +
    scale_color_manual(name = NULL, values = c("SEQ" = "blue", "pi" = "red"), 
                        labels = c("SEQ" = "SEQ",
                                    "pi" = "Default")) +
    scale_y_continuous(labels = function(x) format(round(x, 1), nsmall = 1)) +
    scale_x_continuous(labels = function(x) format(round(x, 1), nsmall = 1)
                        # limits = c(theta_true[3] - c[1], theta_true[3] + c[1])
                      ) +
    theme_minimal()
  plot_s2 <- ggplot() +
    geom_density(aes(x = draws_s2_nppseq, color = "SEQ"), size = 1) +
    geom_density(aes(x = draws_s2_npp, color = "pi"), size = 1) +
    geom_vline(xintercept = theta_true[4], linetype = "dashed", color = "black") +
    labs(x = expression(sigma^2), y = "") +
    scale_color_manual(name = NULL, values = c("SEQ" = "blue", "pi" = "red"), 
                        labels = c("SEQ" = "SEQ",
                                    "pi" = "Default")) +
    scale_y_continuous(labels = function(x) format(round(x, 1), nsmall = 1)) +
    scale_x_continuous(labels = function(x) format(round(x, 1), nsmall = 1),
                        limits = c(0, theta_true[4] + 2)
                      ) +
    theme_minimal() 
  
  out <- plot_beta1 + plot_beta2 + plot_beta3 + plot_s2 +
  plot_layout(ncol = 4, guides = "collect")
  if (leg) {
    out <- out & theme(legend.position = "bottom")
  } else {
    out <- out & theme(legend.position = "none")
  }
  out
}

plot_eta_lm <- function(
  draws, 
  ordered = F,
  post = F,
  title = NULL,
  leg = T
) {
  draws_eta_nppseq <- draws$delta_nppseq[, 1:3]
  draws_eta_npp <- draws$delta_npp[, 1:3]
  if (ordered) {
    draws_eta_nppseq <- draws$delta_onppseq[, 1:3]
    draws_eta_npp <- draws$delta_onpp[, 1:3]
  }
  if (post) {
    y_lab <- expression(paste(pi[A], "(", eta ~ "|" ~ D ~ "," ~ D[0], ")"))
  } else {
    y_lab <- expression(paste(pi[A], "(", eta, ")"))
  }

  plot_eta1 <- ggplot() +
    geom_density(bounds = c(0, 1), aes(x = draws_eta_nppseq[, 1], color = "SEQ"), size = 1) +
    geom_density(bounds = c(0, 1), aes(x = draws_eta_npp[, 1], color = "pi"), size = 1) +
    labs(x = expression(eta[1]), y = y_lab) +
    scale_color_manual(name = NULL, values = c("SEQ" = "blue", "pi" = "red"), 
                        labels = c("SEQ" = "SEQ",
                                    "pi" = "Default")) +
    scale_y_continuous(labels = function(x) format(round(x, 2), nsmall = 2)) +
    scale_x_continuous(labels = function(x) format(round(x, 2), nsmall = 2),
                        limits = c(0,1)) +
    theme_minimal() +
    ggtitle(
      title
    )
  plot_eta2 <- ggplot() +
    geom_density(bounds = c(0, 1), aes(x = draws_eta_nppseq[, 2], color = "SEQ"), size = 1) +
    geom_density(bounds = c(0, 1), aes(x = draws_eta_npp[, 2], color = "pi"), size = 1) +
    labs(x = expression(eta[2]), y = "") +
    scale_color_manual(name = NULL, values = c("SEQ" = "blue", "pi" = "red"), 
                        labels = c("SEQ" = "SEQ",
                                    "pi" = "Default")) +
    scale_y_continuous(labels = function(x) format(round(x, 2), nsmall = 2)) +
    scale_x_continuous(labels = function(x) format(round(x, 2), nsmall = 2),
                        limits = c(0,1)) +
    theme_minimal()
  plot_eta3 <- ggplot() +
    geom_density(bounds = c(0, 1), aes(x = draws_eta_nppseq[, 3], color = "SEQ"), size = 1) +
    geom_density(bounds = c(0, 1), aes(x = draws_eta_npp[, 3], color = "pi"), size = 1) +
    labs(x = expression(eta[3]), y = "") +
    scale_color_manual(name = NULL, values = c("SEQ" = "blue", "pi" = "red"), 
                        labels = c("SEQ" = "SEQ",
                                    "pi" = "Default")) +
    scale_y_continuous(labels = function(x) format(round(x, 2), nsmall = 2)) +
    scale_x_continuous(labels = function(x) format(round(x, 2), nsmall = 2),
                        limits = c(0,1)) +
    theme_minimal()
  
  out <- plot_eta1 + plot_eta2 + plot_eta3 +
  plot_layout(ncol = 3, guides = "collect") 
  if (leg) {
    out <- out & theme(legend.position = "bottom")
  } else {
    out <- out & theme(legend.position = "none")
  }
  out
}

plots_npp_sce <- lapply(1:3, function(i) {
  samples <- sample_prior_post_lm(i, r, sce_list, sim_sces_n0_ge_n) 
  sce_name <- paste0("Scenario ", as.roman(round((i+1)/3)), ".", 
                      as.roman(ifelse(i %% 3 == 0, 3, i %% 3)),
                      ":"
                    )
  title <- bquote(
    .(sce_name)~epsilon == .(df_sce$Eps[i])~", Order:"~.(df_sce$Order[i])~", Discrepancy:"~.(df_sce$Discrepancy[i])
  )
  if (i != 3){
    plot_theta_prior <- plot_theta_lm(samples$theta_prior, title = title, leg = F)
    plot_eta_prior <- plot_eta_lm(samples$eta_prior, title = title, leg = F)
    plot_theta_post <- plot_theta_lm(samples$theta_post, post = T, leg = F)
    plot_eta_post <- plot_eta_lm(samples$eta_post, post = T, leg = F)
  } else {
    plot_theta_prior <- plot_theta_lm(samples$theta_prior, title = title, leg = F)
    plot_eta_prior <- plot_eta_lm(samples$eta_prior, title = title, leg = F)
    plot_theta_post <- plot_theta_lm(samples$theta_post, post = T, leg = T)
    plot_eta_post <- plot_eta_lm(samples$eta_post, post = T, leg = T)
  }
  
  plot_theta <- plot_theta_prior / plot_theta_post
  plot_eta <- plot_eta_prior / plot_eta_post + plot_layout(guides = "collect") 
  list(
    plot_theta = plot_theta,
    plot_eta = plot_eta
    )
})

plots_eta_npp <- lapply(1:3, function(i) plots_npp_sce[[i]]$plot_eta) %>% 
  wrap_plots(ncol = 1, guides = "collect") 
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
