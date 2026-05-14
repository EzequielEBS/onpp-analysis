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
delta_model_lm <- cmdstan_model("code/models/lm/delta_lm.stan")

load("../../simulated-data-regression/data/true_params.RData")
betastar <- beta[[1]]
sgstar <- 1

sce1.1 <- list(
  par = list(
    p = 3,
    a = 2,
    b = 1,
    V0 = diag(3),
    mu0 = rep(0, 3),
    tilde_a = 1,
    tilde_b = 1,
    alpha = rep(1, 4),
    post = 0
  ),
  data = qs_read("results/sim_data/lm/sceI_I.qs2")$data[1:2]
)
sim1.1 <- qs_read("results/samples/lm/sim1_1_lm.qs2")
sim1.2 <- qs_read("results/samples/lm/sim1_2_lm.qs2")
prior_samples <- sim_sce(model = "lm", 15,  1, sce1.1, gamma_model_lm, delta_model_lm)
eta_prior <- prior_samples$delta[[1]]
theta_prior <- prior_samples$theta[[1]]

r <- 3
theta_post_1.1 <- sim1.1$theta[[r]]
theta_post_1.2 <- sim1.2$theta[[r]]
eta_post_1.1 <- sim1.1$delta[[r]]
eta_post_1.2 <- sim1.2$delta[[r]]

mse_1.1 <- lapply(1:4, function(i) {
  theta_mat <- as.matrix(theta_post_1.1[[i]])
  theta_hat <- colMeans(theta_mat)
  true_value <- c(betastar, sgstar)
  mean((theta_hat - true_value)^2)
})
mse_1.2 <- lapply(1:4, function(i) {
  theta_mat <- as.matrix(theta_post_1.2[[i]])
  theta_hat <- colMeans(theta_mat)
  true_value <- c(betastar, sgstar)
  mean((theta_hat - true_value)^2)
})

mse_df <- data.frame(
  prior = c("NPP", "NPP-SEQ", "ONPP", "ONPP-SEQ"),
  mse_sce_1.1 = unlist(mse_1.1),
  mse_sce_1.2 = unlist(mse_1.2)
)
mse_df


plot_theta_lm <- function(draws, theta_true = c(betastar, sgstar), c = c(3, 3), title = NULL) {
  draws_beta_nppseq <- draws$theta_nppseq[, 1:3]
  draws_beta_npp <- draws$theta_npp[, 1:3]
  draws_s2_nppseq <- draws$theta_nppseq[, 4]
  draws_s2_npp <- draws$theta_npp[, 4]

  plot_beta1 <- ggplot() +
    geom_density(aes(x = draws_beta_nppseq[, 1], color = "SEQ"), size = 1) +
    geom_density(aes(x = draws_beta_npp[, 1], color = "pi"), size = 1) +
    labs(x = expression(beta[1]), y = "") +
    geom_vline(xintercept = theta_true[1], linetype = "dashed", color = "black") +
    scale_color_manual(name = NULL, values = c("SEQ" = "blue", "pi" = "red"), 
                        labels = c("SEQ" = "SEQ",
                                    "pi" = "Default")) +
    scale_y_continuous(labels = function(x) format(round(x, 1), nsmall = 1)) +
    scale_x_continuous(labels = function(x) format(round(x, 1), nsmall = 1),
                        limits = c(theta_true[1] - c[1], theta_true[1] + c[1])) +
    theme_minimal() +
    ggtitle(
      title
    )
  plot_beta2 <- ggplot() +
    geom_density(aes(x = draws_beta_nppseq[, 2], color = "SEQ"), size = 1) +
    geom_density(aes(x = draws_beta_npp[, 2], color = "pi"), size = 1) +
    labs(x = expression(beta[2]), y = "") +
    geom_vline(xintercept = theta_true[2], linetype = "dashed", color = "black") +
    scale_color_manual(name = NULL, values = c("SEQ" = "blue", "pi" = "red"), 
                        labels = c("SEQ" = "SEQ",
                                    "pi" = "Default")) +
    scale_y_continuous(labels = function(x) format(round(x, 1), nsmall = 1)) +
    scale_x_continuous(labels = function(x) format(round(x, 1), nsmall = 1),
                        limits = c(theta_true[2] - c[1], theta_true[2] + c[1])) +
    theme_minimal()
  plot_beta3 <- ggplot() +
    geom_density(aes(x = draws_beta_nppseq[, 3], color = "SEQ"), size = 1) +
    geom_density(aes(x = draws_beta_npp[, 3], color = "pi"), size = 1) +
    geom_vline(xintercept = theta_true[3], linetype = "dashed", color = "black") +
    labs(x = expression(beta[3]), y = "") +
    scale_color_manual(name = NULL, values = c("SEQ" = "blue", "pi" = "red"), 
                        labels = c("SEQ" = "SEQ",
                                    "pi" = "Default")) +
    scale_y_continuous(labels = function(x) format(round(x, 1), nsmall = 1)) +
    scale_x_continuous(labels = function(x) format(round(x, 1), nsmall = 1),
                        limits = c(theta_true[3] - c[1], theta_true[3] + c[1])) +
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
                        limits = c(theta_true[4] - c[2], theta_true[4] + c[2])) +
    theme_minimal() 
  
  plot_beta1 + plot_beta2 + plot_beta3 + plot_s2 +
  plot_layout(ncol = 4, guides = "collect") & 
  theme(legend.position = "bottom")
}

plot_eta_lm <- function(draws, title = NULL) {
  draws_eta_nppseq <- draws$delta_nppseq[, 1:3]
  draws_eta_npp <- draws$delta_npp[, 1:3]

  plot_eta1 <- ggplot() +
    geom_density(bounds = c(0, 1), aes(x = draws_eta_nppseq[, 1], color = "SEQ"), size = 1) +
    geom_density(bounds = c(0, 1), aes(x = draws_eta_npp[, 1], color = "pi"), size = 1) +
    labs(x = expression(eta[1]), y = "") +
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
  
  plot_eta1 + plot_eta2 + plot_eta3 +
  plot_layout(ncol = 3, guides = "collect") & 
  theme(legend.position = "bottom")
}

plot_theta_prior <- plot_theta_lm(theta_prior, c = c(4,1), title = "Prior")
plot_eta_prior <- plot_eta_lm(eta_prior, title = "Prior")
plot_theta_post_1.1 <- plot_theta_lm(theta_post_1.1, c = c(0.5,1), title = "Posterior")
plot_eta_post_1.1 <- plot_eta_lm(eta_post_1.1, title = "Posterior")
plot_theta_post_1.2 <- plot_theta_lm(theta_post_1.2, c = c(0.5,1), title = "Posterior")
plot_eta_post_1.2 <- plot_eta_lm(eta_post_1.2, title = "Posterior")

plot_theta_1.1 <- plot_theta_prior / plot_theta_post_1.1 + plot_layout(guides = "collect") & 
  theme(legend.position = "bottom")
plot_eta_1.1 <- plot_eta_prior / plot_eta_post_1.1 + plot_layout(guides = "collect") & 
  theme(legend.position = "bottom")
plot_theta_1.2 <- plot_theta_prior / plot_theta_post_1.2 + plot_layout(guides = "collect") & 
  theme(legend.position = "bottom")
plot_eta_1.2 <- plot_eta_prior / plot_eta_post_1.2 + plot_layout(guides = "collect") & 
  theme(legend.position = "bottom")


ggsave("results/figures/lm/theta_prior_post_npp_sce1_1_lm.png", plot_theta_1.1, width = 12, height = 8)
ggsave("results/figures/lm/eta_prior_post_npp_sce1_1_lm.png", plot_eta_1.1, width = 9, height = 8)
ggsave("results/figures/lm/theta_prior_post_npp_sce1_2_lm.png", plot_theta_1.2, width = 12, height = 8)
ggsave("results/figures/lm/eta_prior_post_npp_sce1_2_lm.png", plot_eta_1.2, width = 9, height = 8)

