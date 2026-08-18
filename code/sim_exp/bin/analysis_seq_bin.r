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
sim1.1 <- sim_sce(model = "bin", n_cores,  n_sim, sce1.1_unif, gamma_model_bin, delta_model_bin)
qs2::qs_save(sim1.1, file = "results/samples/bin/sim1_1_bin_unif.qs2")
sim1.2 <- sim_sce(model = "bin", n_cores,  n_sim, sce1.2_unif, gamma_model_bin, delta_model_bin)
qs2::qs_save(sim1.2, file = "results/samples/bin/sim1_2_bin_unif.qs2")
sim1.3 <- sim_sce(model = "bin", n_cores,  n_sim, sce1.3_unif, gamma_model_bin, delta_model_bin)
qs2::qs_save(sim1.3, file = "results/samples/bin/sim1_3_bin_unif.qs2")

# NOTE: sim1.1/1.2/1.3 already hold exactly what qs_save() just wrote --
# immediately qs_read()-ing them back in was a pure round-trip through disk
# with no effect other than costing the (de)serialization time. Removed. If
# you want to skip re-running the (expensive) sim_sce() calls above on a
# later run, comment them out and uncomment the three qs_read() calls below
# instead of running both unconditionally.
# sim1.1 <- qs2::qs_read("results/samples/bin/sim1_1_bin_unif.qs2")
# sim1.2 <- qs2::qs_read("results/samples/bin/sim1_2_bin_unif.qs2")
# sim1.3 <- qs2::qs_read("results/samples/bin/sim1_3_bin_unif.qs2")

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
  prior_fit <- sample_sce_bin(sce, gamma_model_bin, delta_model_bin)
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


plot_theta_bin <- function(
  draws, 
  ordered = F,
  post = F,
  theta_true = .5, 
  # c = c(3, 3), 
  title = NULL,
  leg = T
) {
  draws_theta_nppseq <- draws$theta_nppseq
  draws_theta_npp <- draws$theta_npp
  if (ordered) {
    draws_theta_nppseq <- draws$theta_onppseq
    draws_theta_npp <- draws$theta_onpp
  }

  if (post) {
    y_lab <- expression(paste(pi[0], "(", theta ~ "|" ~ D ~ "," ~ D[0], ")"))
  } else {
    y_lab <- expression(paste(pi[0], "(", theta, ")"))
  }

  plot_theta <- ggplot() +
    geom_density(aes(x = draws_theta_nppseq, color = "SEQ"), size = 1) +
    geom_density(aes(x = draws_theta_npp, color = "pi"), size = 1) +
    labs(x = expression(theta), y = y_lab) +
    geom_vline(xintercept = theta_true, linetype = "dashed", color = "black") +
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
  if (leg) {
    out <- plot_theta & theme(legend.position = "bottom")
  } else {
    out <- plot_theta & theme(legend.position = "none")
  }
  out
}

plot_eta_bin <- function(
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

r <- 1
plots_npp_sce <- lapply(1:3, function(i) {
  samples <- sample_prior_post_bin(i, r, sce_list, sim_list) 
  sce_name <- paste0("Scenario ", as.roman(round((i+1)/3)), ".", 
                      as.roman(ifelse(i %% 3 == 0, 3, i %% 3)),
                      ":"
                    )
  title <- bquote(
    .(sce_name)~theta[0] == .(df_sce$Theta0[i])~", Order:"~.(df_sce$Order[i])~", Discrepancy:"~.(df_sce$Discrepancy[i])
  )
  if (i != 3){
    plot_theta_prior <- plot_theta_bin(samples$theta_prior, title = title, leg = F)
    plot_eta_prior <- plot_eta_bin(samples$eta_prior, title = title, leg = F)
    plot_theta_post <- plot_theta_bin(samples$theta_post, post = T, leg = F)
    plot_eta_post <- plot_eta_bin(samples$eta_post, post = T, leg = F)
  } else {
    plot_theta_prior <- plot_theta_bin(samples$theta_prior, title = title, leg = F)
    plot_eta_prior <- plot_eta_bin(samples$eta_prior, title = title, leg = F)
    plot_theta_post <- plot_theta_bin(samples$theta_post, post = T, leg = T)
    plot_eta_post <- plot_eta_bin(samples$eta_post, post = T, leg = T)
  }
  
  plot_theta <- plot_theta_prior + plot_theta_post
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

ggsave("results/figures/bin/theta_prior_post_npp.pdf", plots_theta_npp, width = 9, height = 13.5)
ggsave("results/figures/bin/eta_prior_post_npp.pdf", plots_eta_npp, width = 9, height = 13.5)
