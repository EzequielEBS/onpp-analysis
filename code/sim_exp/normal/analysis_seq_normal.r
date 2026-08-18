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
gamma_model_normal <- cmdstan_model("code/models/normal/normal_fixed_var_onpp.stan")
delta_model_normal <- cmdstan_model("code/models/normal/normal_fixed_var_npp.stan")


# define scenarios
## compatible order
### compatible order and high congruence
sce1.1 <- list(n0 = c(100, 100, 100),
              n = 100,
              sigma2 = 1,
              mu0 = 0,
              tau02 = 1,
              al = 1/2,
              bl = 1/2,
              theta0 = c(4.6, 4.8, 5),
              theta = 5,
              alpha = rep(1/4, 4),
              post = 0
)
### compatible order and small congruence
sce1.2 <- list(n0 = c(100, 100, 100),
              n = 100,
              sigma2 = 1,
              mu0 = 0,
              tau02 = 1,
              al = 1/2,
              bl = 1/2,
              theta0 = c(3.6, 3.8, 4),
              theta = 5,
              alpha = rep(1/4, 4),
              post = 0
)
### compatible order and no congruence 
sce1.3 <- list(n0 = c(100, 100, 100),
              n = 100,
              sigma2 = 1,
              mu0 = 0,
              tau02 = 1,
              al = 1/2,
              bl = 1/2,
              theta0 = c(2.6, 2.8, 3),
              theta = 5,
              alpha = rep(1/4, 4),
              post = 0
)


# load sample
sim1.1 <- qs_read("results/samples/normal/sim1_1_normal_fixed_var.qs2")
sim1.2 <- qs_read("results/samples/normal/sim1_2_normal_fixed_var.qs2")
sim1.3 <- qs_read("results/samples/normal/sim1_3_normal_fixed_var.qs2")

sce_list <- list(sce1.1, sce1.2, sce1.3)
sim_list <- list(sim1.1, sim1.2, sim1.3)

df_sce <- data.frame(
  # Regime = c(rep("n[0] > n", 12), rep("n[0] < n", 12)),
  Theta0 = c(
    "(4.6, 4.8, 5.0)",
    "(3.6, 3.8, 4.0)",
    "(2.6, 2.8, 3.0)",
    "(4.6, 5.0, 4.8)",
    "(3.6, 4.0, 3.8)",
    "(2.6, 3.0, 2.8)",
    "(5.0, 4.8, 4.6)",
    "(4.0, 3.8, 3.6)",
    "(3.0, 2.8, 2.6)",
    "(5.0, 5.0, 5.0)",
    "(4.0, 4.0, 4.0)",
    "(3.0, 3.0, 3.0)",
    "(3.0, 7.0, 5.0)",
    "(5.0, 7.0, 3.0)"
  ),
  Order       = c(rep("Compatible",    3),
                  rep("Semi-compat.",  3),
                  rep("Incompatible",  3),
                  rep("Neutral",       3),
                  rep("Cancellation effect", 2)
                  ),
  Discrepancy = c(rep(c("None","Small","Large"), 4), rep("Large", 2))
                )

sample_prior_post_normal <- function(i, r, sce_list, sim_list) {
  sce <- sce_list[[i]]
  sim <- sim_list[[i]]
  # This only ever needs ONE prior draw (only prior_samples$delta[[1]] /
  # $theta[[1]] were used), but used to get it via sim_sce(model=
  # "normal_fixed_var", 2, 2, ...) -- which spins up a 2-worker PSOCK
  # cluster, reloads cmdstanr/dplyr/posterior on both workers, fits 2 full
  # replicates, and then discards the second one. sample_sce_normal_fixed_var()
  # is the same per-replicate fit sim_sce() calls internally; calling it
  # directly here (sce$post is already 0, so this draws from the prior)
  # gives the same single draw without any cluster setup/teardown overhead.
  # Same fix as analysis_seq_bin.r's sample_prior_post_bin() and
  # analysis_seq_lm.r's sample_prior_post_lm().
  prior_fit <- sample_sce_normal_fixed_var(sce, gamma_model_normal, delta_model_normal)
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


plot_theta_normal <- function(
  draws, 
  ordered = F,
  post = F,
  theta_true = 5, 
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

plot_eta_normal <- function(
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
  samples <- sample_prior_post_normal(i, r, sce_list, sim_list) 
  sce_name <- paste0("Scenario ", as.roman(round((i+1)/3)), ".", 
                      as.roman(ifelse(i %% 3 == 0, 3, i %% 3)),
                      ":"
                    )
  title <- bquote(
    .(sce_name)~theta[0] == .(df_sce$Theta0[i])~", Order:"~.(df_sce$Order[i])~", Discrepancy:"~.(df_sce$Discrepancy[i])
  )
  if (i != 3){
    plot_theta_prior <- plot_theta_normal(samples$theta_prior, title = title, leg = F)
    plot_eta_prior <- plot_eta_normal(samples$eta_prior, title = title, leg = F)
    plot_theta_post <- plot_theta_normal(samples$theta_post, post = T, leg = F)
    plot_eta_post <- plot_eta_normal(samples$eta_post, post = T, leg = F)
  } else {
    plot_theta_prior <- plot_theta_normal(samples$theta_prior, title = title, leg = F)
    plot_eta_prior <- plot_eta_normal(samples$eta_prior, title = title, leg = F)
    plot_theta_post <- plot_theta_normal(samples$theta_post, post = T,leg = T)
    plot_eta_post <- plot_eta_normal(samples$eta_post, post = T, leg = T)
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

ggsave("results/figures/normal/theta_prior_post_npp.pdf", plots_theta_npp, width = 9, height = 13.5)
ggsave("results/figures/normal/eta_prior_post_npp.pdf", plots_eta_npp, width = 9, height = 13.5)
