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

# Compile the model
gamma_model_bin <- cmdstan_model("code/models/bin/gamma_bin.stan")
delta_model_bin <- cmdstan_model("code/models/bin/eta_bin.stan")

# define scenarios
#
# NOTE: these use a = b = 1 (a flat/uniform gamma prior), unlike the
# same-numbered scenarios in simulations_bin.R's scenario table (a = b =
# 1/2). They are named sceX.Y_unif / simX.Y_unif here -- matching the
# "_unif" suffix already used for their output files -- rather than reusing
# the bare sce1.1/sim1.1 names, so they can't be confused with (or
# accidentally shadow) the differently-parameterized default-prior versions
# defined in simulations_bin.R.
## compatible order
### compatible order and high congruence
sce1.1_unif <- list(n0 = c(100, 100, 100),
              n = 100,
              a = 1,
              b = 1,
              al = 1/2,
              bl = 1/2,
              theta0 = c(0.3, 0.4, 0.5),
              theta = 0.5,
              alpha = rep(1/4, 4),
              post = 1
)
### compatible order and small congruence
sce1.2_unif <- list(n0 = c(100, 100, 100),
             n = 100,
             a = 1,
             b = 1,
             al = 1/2,
             bl = 1/2,
             theta0 = c(0.2, 0.3, 0.4),
             theta = 0.5,
             alpha = rep(1/4, 4),
              post = 1
)
### compatible order and no congruence
sce1.3_unif <- list(n0 = c(100, 100, 100),
             n = 100,
             a = 1,
             b = 1,
             al = 1/2,
             bl = 1/2,
             theta0 = c(0.05, 0.15, 0.25),
             theta = 0.5,
             alpha = rep(1/4, 4),
              post = 1
)


# load sample
n_cores <- 15
n_sim <- 200
# sim1.1 <- sim_sce(model = "bin", n_cores,  n_sim, sce1.1_unif, gamma_model_bin, delta_model_bin)
# qs2::qs_save(sim1.1, file = "results/samples/bin/sim1_1_bin_unif.qs2")
# sim1.2 <- sim_sce(model = "bin", n_cores,  n_sim, sce1.2_unif, gamma_model_bin, delta_model_bin)
# qs2::qs_save(sim1.2, file = "results/samples/bin/sim1_2_bin_unif.qs2")
# sim1.3 <- sim_sce(model = "bin", n_cores,  n_sim, sce1.3_unif, gamma_model_bin, delta_model_bin)
# qs2::qs_save(sim1.3, file = "results/samples/bin/sim1_3_bin_unif.qs2")

# NOTE: sim1.1/1.2/1.3 already hold exactly what qs_save() just wrote --
# immediately qs_read()-ing them back in was a pure round-trip through disk
# with no effect other than costing the (de)serialization time. Removed. If
# you want to skip re-running the (expensive) sim_sce() calls above on a
# later run, comment them out and uncomment the three qs_read() calls below
# instead of running both unconditionally.
sim1.1 <- qs2::qs_read("results/samples/bin/sim1_1_bin_unif.qs2")
sim1.2 <- qs2::qs_read("results/samples/bin/sim1_2_bin_unif.qs2")
sim1.3 <- qs2::qs_read("results/samples/bin/sim1_3_bin_unif.qs2")

sce1.1_unif$post <- 0
sce1.2_unif$post <- 0
sce1.3_unif$post <- 0
sce_list <- list(sce1.1_unif, sce1.2_unif, sce1.3_unif)
sim_list <- list(sim1.1, sim1.2, sim1.3)

df_sce <- data.frame(
  # Regime = c(rep("n[0] > n", 12), rep("n[0] < n", 12)),
  Theta0 = c(
    "(0.30, 0.40, 0.50)",
    "(0.20, 0.30, 0.40)",
    "(0.05, 0.15, 0.25)",
    "(0.30, 0.50, 0.40)",
    "(0.20, 0.40, 0.30)",
    "(0.05, 0.25, 0.15)",
    "(0.50, 0.40, 0.30)",
    "(0.40, 0.30, 0.20)",
    "(0.25, 0.15, 0.05)",
    "(0.50, 0.50, 0.50)",
    "(0.30, 0.30, 0.30)",
    "(0.10, 0.10, 0.10)",
    "(0.10, 0.90, 0.50)",
    "(0.50, 0.90, 0.10)"
  ),
  Order       = c(rep("Compatible",    3),
                  rep("Semi-compat.",  3),
                  rep("Incompatible",  3),
                  rep("Neutral",       3),
                  rep("Cancellation effect", 2)
                  ),
  Discrepancy = c(rep(c("None","Small","Large"), 4), rep("Large", 2))
                )

sample_prior_post_bin <- function(i, r, sce_list, sim_list) {
  sce <- sce_list[[i]]
  sim <- sim_list[[i]]
  # This only ever needs ONE prior draw (only prior_samples$delta[[1]] /
  # $theta[[1]] were used), but used to get it via sim_sce(model="bin", 2,
  # 2, ...) -- which spins up a 2-worker PSOCK cluster, reloads cmdstanr/
  # dplyr/posterior on both workers, fits 2 full replicates, and then
  # discards the second one. sample_sce_bin() is the same per-replicate fit
  # sim_sce() calls internally; calling it directly here (sce$post is
  # already forced to 0 above, so this draws from the prior) gives the same
  # single draw without any cluster setup/teardown overhead.
  #
  # Cached to disk (qs2) so re-running this script -- e.g. just to tweak
  # figure aesthetics -- doesn't re-fit the Stan prior model every time.
  # Delete the cache file to force a re-sample (e.g. after changing a
  # scenario's parameters in sce_list above).
  cache_file <- sprintf("results/samples/bin/prior_cache_sce%d.qs2", i)
  if (file.exists(cache_file)) {
    prior_fit <- qs2::qs_read(cache_file)
  } else {
    prior_fit <- sample_sce_bin(sce, gamma_model_bin, delta_model_bin)
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


# plot_theta_bin()/plot_eta_bin() used to be defined here (and, byte-for-
# byte identically apart from the family name, in analysis_seq_normal.r and
# analysis_seq_lm.r). They're now plot_theta_prior_post()/
# plot_eta_prior_post() in aux_fun_sim.R -- one shared definition instead
# of three copies that could silently drift apart -- see the comment above
# those functions for the full rationale.

r <- 1
plots_npp_sce <- lapply(1:3, function(i) {
  samples <- sample_prior_post_bin(i, r, sce_list, sim_list)
  # No "Scenario ..." label -- theta0/Order/Discrepancy already identify
  # the scenario.
  title <- bquote(
    theta[0] == .(df_sce$Theta0[i])~", Order:"~.(df_sce$Order[i])~", Discrepancy:"~.(df_sce$Discrepancy[i])
  )
  if (i != 3){
    plot_theta_prior <- plot_theta_prior_post(samples$theta_prior, title = title, leg = F,
      theta_true = 0.5, low_band = 0.5, upp_band = 0.5)
    plot_eta_prior <- plot_eta_prior_post(samples$eta_prior, title = title, leg = F)
    plot_theta_post <- plot_theta_prior_post(samples$theta_post, post = T, leg = F,
      theta_true = 0.5, low_band = 0.2, upp_band = 0.2)
    plot_eta_post <- plot_eta_prior_post(samples$eta_post, post = T, leg = F)
  } else {
    plot_theta_prior <- plot_theta_prior_post(samples$theta_prior, title = title, leg = F,
      theta_true = 0.5, low_band = 0.5, upp_band = 0.5)
    plot_eta_prior <- plot_eta_prior_post(samples$eta_prior, title = title, leg = F)
    plot_theta_post <- plot_theta_prior_post(samples$theta_post, post = T, leg = T,
      theta_true = 0.5, low_band = 0.2, upp_band = 0.2)
    plot_eta_post <- plot_eta_prior_post(samples$eta_post, post = T, leg = T)
  }
  
  plot_theta <- plot_theta_prior + plot_theta_post
  # No plot_layout(guides = "collect") here -- see the note below on
  # plots_eta_npp for why relying on collect (rather than the leg flag
  # alone) was producing duplicate/misplaced legends.
  plot_eta <- plot_eta_prior / plot_eta_post
  list(
    plot_theta = plot_theta,
    plot_eta = plot_eta
    )
})

# NOTE on legends: this composition deliberately avoids
# plot_layout(guides = "collect")/wrap_plots(..., guides = "collect")
# entirely for the eta figure -- exactly like plots_theta_npp below,
# which never used collect and has always rendered with a single legend.
# Two things were tried and both produced duplicate/misplaced legends:
#  1. guides = "collect" at both the per-scenario level (plot_eta_prior /
#     plot_eta_post) and again at this outer wrap_plots() -- once
#     scenario 3's own local collect had already resolved/materialized
#     its legend, this outer collect couldn't merge it with the (empty)
#     guide sets of scenarios 1-2, and just left it where the inner
#     collect put it, positioned by patchwork's default ("right") since
#     no theme(legend.position=...) was set at this true top level.
#  2. Adding `& theme(legend.position = "bottom")` here to fix that
#     positioning: `&` recurses into every nested sub-plot, including the
#     ones with legend.position = "none" from leg = F at the call sites
#     above -- forcing position = "bottom" on those rebuilt a guide for
#     each of them too (3 separate legends instead of 1).
# Instead: only scenario 3's post plot ever sets leg = T, so it's the
# only one that builds a guide at all; with no collect anywhere, that one
# guide renders right where plot_eta_prior_post() already puts it (its
# own theme(legend.position = "bottom")), and since it's the last block,
# that's the bottom of the whole figure.
plots_eta_npp <- lapply(1:3, function(i) plots_npp_sce[[i]]$plot_eta) %>%
  wrap_plots(ncol = 1)
plots_theta_npp <- lapply(1:3, function(i) plots_npp_sce[[i]]$plot_theta) %>%
  wrap_plots(ncol = 1)

ggsave("results/figures/bin/theta_prior_post_npp.pdf", plots_theta_npp, width = 9, height = 13.5)
ggsave("results/figures/bin/eta_prior_post_npp.pdf", plots_eta_npp, width = 9, height = 13.5)
