# Load data
load("results/sim_data/perf_sce_lm_data.RData")
load("results/sim_data/inc_sce_lm_data.RData")

# Load libraries
library(cmdstanr)
library(posterior)
library(tidyverse)
library(dplyr)
library(mvtnorm)
library(ggplot2)
library(parallel)
library(MCMCpack)
library(ggh4x)
library(rstan)
library(shinystan)
library(hdbayes)
library(parallel)

# Compile the model
gamma_model <- cmdstan_model("code/models/gamma_lm.stan")
delta_model <- cmdstan_model("code/models/delta_lm.stan")

fit_lm <- function(data, models, prior_parameters) {
  # Get prior parameters
  a <- prior_parameters$a
  b <- prior_parameters$b
  tilde_a <- prior_parameters$tilde_a
  tilde_b <- prior_parameters$tilde_b
  alpha <- prior_parameters$alpha
  mu0 <- prior_parameters$mu0
  V0 <- prior_parameters$V0
  
  # Prepare data for Stan
  stan_data_npp <- list(
    K = nrow(data$dims_X0),
    n0 =data$dims_X0[,1],
    n = nrow(data$X),
    p = ncol(data$X),
    len_X0 = length(data$X0_flat),
    len_y0 = length(data$y0_flat),
    X0 = data$X0_flat,
    startid_X0 = data$startid_X0,
    y0 = data$y0_flat,
    startid_y0 = data$startid_y0,
    X = data$X,
    y = as.vector(data$y),
    a = a,
    b = b,
    V0 = V0,
    mu0 = mu0,
    tilde_a = tilde_a,
    tilde_b = tilde_b,
    post = 1,
    seq = 0
  )
  stan_data_nppseq <- list(
    K = nrow(data$dims_X0),
    n0 =data$dims_X0[,1],
    n = nrow(data$X),
    p = ncol(data$X),
    len_X0 = length(data$X0_flat),
    len_y0 = length(data$y0_flat),
    X0 = data$X0_flat,
    startid_X0 = data$startid_X0,
    y0 = data$y0_flat,
    startid_y0 = data$startid_y0,
    X = data$X,
    y = as.vector(data$y),
    a = a,
    b = b,
    V0 = V0,
    mu0 = mu0,
    tilde_a = tilde_a,
    tilde_b = tilde_b,
    post = 1,
    seq = 1
  )
  stan_data_onpp <- list(
    K = nrow(data$dims_X0),
    n0 =data$dims_X0[,1],
    n = nrow(data$X),
    p = ncol(data$X),
    len_X0 = length(data$X0_flat),
    len_y0 = length(data$y0_flat),
    X0 = data$X0_flat,
    startid_X0 = data$startid_X0,
    y0 = data$y0_flat,
    startid_y0 = data$startid_y0,
    X = data$X,
    y = as.vector(data$y),
    a = a,
    b = b,
    V0 = V0,
    mu0 = mu0,
    alpha = alpha,
    post = 1,
    seq = 0
  )
  stan_data_onppseq <- list(
    K = nrow(data$dims_X0),
    n0 =data$dims_X0[,1],
    n = nrow(data$X),
    p = ncol(data$X),
    len_X0 = length(data$X0_flat),
    len_y0 = length(data$y0_flat),
    X0 = data$X0_flat,
    startid_X0 = data$startid_X0,
    y0 = data$y0_flat,
    startid_y0 = data$startid_y0,
    X = data$X,
    y = as.vector(data$y),
    a = a,
    b = b,
    V0 = V0,
    mu0 = mu0,
    alpha = alpha,
    post = 1,
    seq = 1
  )
  
  gamma_model <- models$gamma_model
  delta_model <- models$delta_model
  
  # Fit the model
  # divergences_onpp <- 1
  # while (divergences_onpp > 0) {
  #   sample_delta_onpp <- model$sample(
  #     data = stan_data_onpp,
  #     chains = 4,
  #     parallel_chains = 4,
  #     iter_warmup = 3000,
  #     iter_sampling = 2000,
  #     refresh = 0,
  #     adapt_delta = 0.99,
  #     step_size = 0.03
  #   )
  #   
  #   # Check for divergences
  #   diagnostics_onpp <- sample_delta_onpp$sampler_diagnostics()
  #   divergences_onpp <- sum((diagnostics_onpp[, , "divergent__"] %>% as_draws_df())$divergent__)
  # }
  sample_delta_npp <- delta_model$sample(
    data = stan_data_npp,
    chains = 4,
    parallel_chains = 4,
    iter_warmup = 4000,
    iter_sampling = 2000,
    refresh = 0,
    adapt_delta = 0.99,
    step_size = 0.03
  )
  sample_delta_nppseq <- delta_model$sample(
    data = stan_data_nppseq,
    chains = 4,
    parallel_chains = 4,
    iter_warmup = 4000,
    iter_sampling = 2000,
    refresh = 0,
    adapt_delta = 0.99,
    step_size = 0.03
  )
  sample_delta_onpp <- gamma_model$sample(
    data = stan_data_onpp,
    chains = 4,
    parallel_chains = 4,
    iter_warmup = 4000,
    iter_sampling = 2000,
    refresh = 0,
    adapt_delta = 0.99,
    step_size = 0.03
  )
  
  sample_delta_onppseq <- gamma_model$sample(
    data = stan_data_onppseq,
    chains = 4,
    parallel_chains = 4,
    iter_warmup = 4000,
    iter_sampling = 2000,
    refresh = 0,
    adapt_delta = 0.99,
    step_size = 0.03
  )
  
  draws_delta_npp <- sample_delta_npp$draws(variables = "delta") %>% 
    as_draws_matrix()
  draws_delta_nppseq <- sample_delta_nppseq$draws(variables = "delta") %>%
    as_draws_matrix()
  draws_delta_onpp <- sample_delta_onpp$draws(variables = "delta") %>% 
    as_draws_matrix()
  draws_delta_onppseq <- sample_delta_onppseq$draws(variables = "delta") %>%
    as_draws_matrix()
  draws_beta_npp <- sample_delta_npp$draws(variables = "beta") %>% 
    as_draws_matrix()
  draws_beta_nppseq <- sample_delta_nppseq$draws(variables = "beta") %>%
    as_draws_matrix()
  draws_beta_onpp <- sample_delta_onpp$draws(variables = "beta") %>% 
    as_draws_matrix()
  draws_beta_onppseq <- sample_delta_onppseq$draws(variables = "beta") %>%
    as_draws_matrix()
  draws_sigma_npp <- sample_delta_npp$draws(variables = "sigma") %>%
    as_draws_matrix()
  draws_sigma_nppseq <- sample_delta_nppseq$draws(variables = "sigma") %>%
    as_draws_matrix()
  draws_sigma_onpp <- sample_delta_onpp$draws(variables = "sigma") %>%
    as_draws_matrix()
  draws_sigma_onppseq <- sample_delta_onppseq$draws(variables = "sigma") %>%
    as_draws_matrix()
  return(list(
    draws_delta_npp = draws_delta_npp,
    draws_delta_nppseq = draws_delta_nppseq,
    draws_delta_onpp = draws_delta_onpp,
    draws_delta_onppseq = draws_delta_onppseq,
    draws_beta_npp = draws_beta_npp,
    draws_beta_nppseq = draws_beta_nppseq,
    draws_beta_onpp = draws_beta_onpp,
    draws_beta_onppseq = draws_beta_onppseq,
    draws_sigma_npp = draws_sigma_npp,
    draws_sigma_nppseq = draws_sigma_nppseq,
    draws_sigma_onpp = draws_sigma_onpp,
    draws_sigma_onppseq = draws_sigma_onppseq
  ))
}


# define prior parameters
K <- length(perf_sce_data[[1]]$X0)
p <- ncol(perf_sce_data[[1]]$X)
perf_sce_hyper <- list(
  a = 2,
  b = 1,
  tilde_a = 1/2,
  tilde_b = 1/2,
  alpha = rep(1/(K+1), K+1),
  mu0 = rep(0, p),
  V0 = diag(p)
)

# Fitting datasets in parallel
sample_perf_sce <- mclapply(perf_sce_data, function(data) {
  fit_lm(data, 
         list(gamma_model = gamma_model,
              delta_model = delta_model),
         perf_sce_hyper)
})
sample_inc_sce <- mclapply(inc_sce_data, function(data) {
  fit_lm(data, 
         list(gamma_model = gamma_model,
              delta_model = delta_model),
         perf_sce_hyper)
})

# Extract draws from the samples
get_hat_par <- function(samples) {
  hat_delta_npp <- do.call(rbind, 
                         lapply(samples, 
                                function(x) colMeans(x$draws_delta_npp)))
  hat_delta_nppseq <- do.call(rbind,
                            lapply(samples, 
                                   function(x) colMeans(x$draws_delta_nppseq)))
  hat_delta_onpp <- do.call(rbind,
                          lapply(samples, 
                                 function(x) colMeans(x$draws_delta_onpp)))
  hat_delta_onppseq <- do.call(rbind,
                             lapply(samples, 
                                    function(x) colMeans(x$draws_delta_onppseq)))
  hat_beta_npp <- do.call(rbind,
                        lapply(samples, 
                               function(x) colMeans(x$draws_beta_npp)))
  hat_beta_nppseq <- do.call(rbind,
                        lapply(samples, 
                               function(x) colMeans(x$draws_beta_nppseq)))
  hat_beta_onpp <- do.call(rbind,
                        lapply(samples, 
                               function(x) colMeans(x$draws_beta_onpp)))
  hat_beta_onppseq <- do.call(rbind,
                        lapply(samples, 
                               function(x) colMeans(x$draws_beta_onppseq)))
  hat_sigma_npp <- do.call(rbind,
                        lapply(samples, 
                               function(x) colMeans(x$draws_sigma_npp)))
  hat_sigma_nppseq <- do.call(rbind,
                        lapply(samples, 
                               function(x) colMeans(x$draws_sigma_nppseq)))
  hat_sigma_onpp <- do.call(rbind,
                        lapply(samples, 
                               function(x) colMeans(x$draws_sigma_onpp)))
  hat_sigma_onppseq <- do.call(rbind,
                        lapply(samples, 
                               function(x) colMeans(x$draws_sigma_onppseq)))
  return(list(hat_delta_npp = hat_delta_npp,
              hat_delta_nppseq = hat_delta_nppseq,
              hat_delta_onpp = hat_delta_onpp,
              hat_delta_onppseq = hat_delta_onppseq,
              hat_beta_npp = hat_beta_npp,
              hat_beta_nppseq = hat_beta_nppseq,
              hat_beta_onpp = hat_beta_onpp,
              hat_beta_onppseq = hat_beta_onppseq,
              hat_sigma_npp = hat_sigma_npp,
              hat_sigma_nppseq = hat_sigma_nppseq,
              hat_sigma_onpp = hat_sigma_onpp,
              hat_sigma_onppseq = hat_sigma_onppseq))
}

# Plot stuff for perfect scenario

hat_par_perf_sce <- get_hat_par(sample_perf_sce)

hat_delta1_perf_sce <- cbind(hat_par_perf_sce$hat_delta_npp[,1], 
                              hat_par_perf_sce$hat_delta_nppseq[,1],
                              hat_par_perf_sce$hat_delta_onpp[,1],
                              hat_par_perf_sce$hat_delta_onppseq[,1])
hat_delta2_perf_sce <- cbind(hat_par_perf_sce$hat_delta_npp[,2],
                              hat_par_perf_sce$hat_delta_nppseq[,2],
                              hat_par_perf_sce$hat_delta_onpp[,2],
                              hat_par_perf_sce$hat_delta_onppseq[,2])
hat_delta3_perf_sce <- cbind(hat_par_perf_sce$hat_delta_npp[,3],
                              hat_par_perf_sce$hat_delta_nppseq[,3],
                              hat_par_perf_sce$hat_delta_onpp[,3],
                              hat_par_perf_sce$hat_delta_onppseq[,3])

hat_beta1_perf_sce <- cbind(hat_par_perf_sce$hat_beta_npp[,1],
                             hat_par_perf_sce$hat_beta_nppseq[,1],
                             hat_par_perf_sce$hat_beta_onpp[,1],
                             hat_par_perf_sce$hat_beta_onppseq[,1])
hat_beta2_perf_sce <- cbind(hat_par_perf_sce$hat_beta_npp[,2],
                             hat_par_perf_sce$hat_beta_nppseq[,2],
                             hat_par_perf_sce$hat_beta_onpp[,2],
                             hat_par_perf_sce$hat_beta_onppseq[,2])
hat_beta3_perf_sce <- cbind(hat_par_perf_sce$hat_beta_npp[,3],
                             hat_par_perf_sce$hat_beta_nppseq[,3],
                             hat_par_perf_sce$hat_beta_onpp[,3],
                             hat_par_perf_sce$hat_beta_onppseq[,3])

hat_sigma_perf_sce <- cbind(hat_par_perf_sce$hat_sigma_npp[,1],
                             hat_par_perf_sce$hat_sigma_nppseq[,1],
                             hat_par_perf_sce$hat_sigma_onpp[,1],
                             hat_par_perf_sce$hat_sigma_onppseq[,1])

# Plotting


box_delta1_perf_sce <- plot_boxplot(hat_delta1_perf_sce, 
                                     expression(hat(delta)[1]))
box_delta2_perf_sce <- plot_boxplot(hat_delta2_perf_sce,
                                     expression(hat(delta)[2]))
box_delta3_perf_sce <- plot_boxplot(hat_delta3_perf_sce,
                                     expression(hat(delta)[3]))

box_beta1_perf_sce <- plot_boxplot(hat_beta1_perf_sce, 
                                      expression(hat(beta)[1])) + 
  geom_hline(yintercept = perf_sce_data[[1]]$beta[1],
             color = "red", linetype = "dashed", size = 1)
box_beta2_perf_sce <- plot_boxplot(hat_beta2_perf_sce, 
                                      expression(hat(beta)[2])) +
  geom_hline(yintercept = perf_sce_data[[1]]$beta[2],
             color = "red", linetype = "dashed", size = 1)
box_beta3_perf_sce <- plot_boxplot(hat_beta3_perf_sce, 
                                      expression(hat(beta)[3])) +
  geom_hline(yintercept = perf_sce_data[[1]]$beta[3],
             color = "red", linetype = "dashed", size = 1)

box_sigma_perf_sce <- plot_boxplot(hat_sigma_perf_sce, 
                                      expression(hat(sigma^2))) +
  geom_hline(yintercept = perf_sce_data[[1]]$sg,
             color = "red", linetype = "dashed", size = 1)


# save plots
ggsave("results/figures/lm_sim/perf_sce/perf_sce_delta1.pdf", 
       plot = box_delta1_perf_sce, 
       width = 8, height = 6)
ggsave("results/figures/lm_sim/perf_sce/perf_sce_delta2.pdf",
       plot = box_delta2_perf_sce, 
       width = 8, height = 6)
ggsave("results/figures/lm_sim/perf_sce/perf_sce_delta3.pdf",
       plot = box_delta3_perf_sce, 
       width = 8, height = 6)
ggsave("results/figures/lm_sim/perf_sce/perf_sce_beta1.pdf",
       plot = box_beta1_perf_sce, 
       width = 8, height = 6)
ggsave("results/figures/lm_sim/perf_sce/perf_sce_beta2.pdf",
       plot = box_beta2_perf_sce, 
       width = 8, height = 6)
ggsave("results/figures/lm_sim/perf_sce/perf_sce_beta3.pdf",
       plot = box_beta3_perf_sce, 
       width = 8, height = 6)
ggsave("results/figures/lm_sim/perf_sce/perf_sce_sigma.pdf",
       plot = box_sigma_perf_sce, 
       width = 8, height = 6)

# -------------------------------------------------------------------------------
# Plotting stuff for incompatible scenario

hat_par_inc_sce <- get_hat_par(sample_inc_sce)
hat_delta1_inc_sce <- cbind(hat_par_inc_sce$hat_delta_npp[,1], 
                             hat_par_inc_sce$hat_delta_nppseq[,1],
                             hat_par_inc_sce$hat_delta_onpp[,1],
                             hat_par_inc_sce$hat_delta_onppseq[,1])
hat_delta2_inc_sce <- cbind(hat_par_inc_sce$hat_delta_npp[,2],
                             hat_par_inc_sce$hat_delta_nppseq[,2],
                             hat_par_inc_sce$hat_delta_onpp[,2],
                             hat_par_inc_sce$hat_delta_onppseq[,2])
hat_delta3_inc_sce <- cbind(hat_par_inc_sce$hat_delta_npp[,3],
                             hat_par_inc_sce$hat_delta_nppseq[,3],
                             hat_par_inc_sce$hat_delta_onpp[,3],
                             hat_par_inc_sce$hat_delta_onppseq[,3])
hat_beta1_inc_sce <- cbind(hat_par_inc_sce$hat_beta_npp[,1],
                            hat_par_inc_sce$hat_beta_nppseq[,1],
                            hat_par_inc_sce$hat_beta_onpp[,1],
                            hat_par_inc_sce$hat_beta_onppseq[,1])
hat_beta2_inc_sce <- cbind(hat_par_inc_sce$hat_beta_npp[,2],
                            hat_par_inc_sce$hat_beta_nppseq[,2],
                            hat_par_inc_sce$hat_beta_onpp[,2],
                            hat_par_inc_sce$hat_beta_onppseq[,2])
hat_beta3_inc_sce <- cbind(hat_par_inc_sce$hat_beta_npp[,3],
                            hat_par_inc_sce$hat_beta_nppseq[,3],
                            hat_par_inc_sce$hat_beta_onpp[,3],
                            hat_par_inc_sce$hat_beta_onppseq[,3])
hat_sigma_inc_sce <- cbind(hat_par_inc_sce$hat_sigma_npp[,1],
                            hat_par_inc_sce$hat_sigma_nppseq[,1],
                            hat_par_inc_sce$hat_sigma_onpp[,1],
                            hat_par_inc_sce$hat_sigma_onppseq[,1])
# Plotting
box_delta1_inc_sce <- plot_boxplot(hat_delta1_inc_sce, 
                                    expression(hat(delta)[1]))
box_delta2_inc_sce <- plot_boxplot(hat_delta2_inc_sce,
                                    expression(hat(delta)[2]))
box_delta3_inc_sce <- plot_boxplot(hat_delta3_inc_sce,
                                    expression(hat(delta)[3]))
box_beta1_inc_sce <- plot_boxplot(hat_beta1_inc_sce,
                                   expression(hat(beta)[1])) + 
  geom_hline(yintercept = inc_sce_data[[1]]$beta[1],
             color = "red", linetype = "dashed", size = 1)
box_beta2_inc_sce <- plot_boxplot(hat_beta2_inc_sce,
                                   expression(hat(beta)[2])) +
  geom_hline(yintercept = inc_sce_data[[1]]$beta[2],
             color = "red", linetype = "dashed", size = 1)
box_beta3_inc_sce <- plot_boxplot(hat_beta3_inc_sce,
                                   expression(hat(beta)[3])) +
  geom_hline(yintercept = inc_sce_data[[1]]$beta[3],
             color = "red", linetype = "dashed", size = 1)
box_sigma_inc_sce <- plot_boxplot(hat_sigma_inc_sce,
                                   expression(hat(sigma^2))) +
  geom_hline(yintercept = inc_sce_data[[1]]$sg,
             color = "red", linetype = "dashed", size = 1)
# save plots
ggsave("results/figures/lm_sim/inc_sce/inc_sce_delta1.pdf",
       plot = box_delta1_inc_sce, 
       width = 8, height = 6)
ggsave("results/figures/lm_sim/inc_sce/inc_sce_delta2.pdf",
       plot = box_delta2_inc_sce, 
       width = 8, height = 6)
ggsave("results/figures/lm_sim/inc_sce/inc_sce_delta3.pdf",
       plot = box_delta3_inc_sce, 
       width = 8, height = 6)
ggsave("results/figures/lm_sim/inc_sce/inc_sce_beta1.pdf",
       plot = box_beta1_inc_sce, 
       width = 8, height = 6)
ggsave("results/figures/lm_sim/inc_sce/inc_sce_beta2.pdf",
       plot = box_beta2_inc_sce, 
       width = 8, height = 6)
ggsave("results/figures/lm_sim/inc_sce/inc_sce_beta3.pdf",
       plot = box_beta3_inc_sce, 
       width = 8, height = 6)
ggsave("results/figures/lm_sim/inc_sce/inc_sce_sigma.pdf",
       plot = box_sigma_inc_sce, 
       width = 8, height = 6)




#-------------------------------------------------------------------------------
# diagostics
#-------------------------------------------------------------------------------

K <- length(sce1_data[[1]]$X0)
p <- ncol(sce1_data[[1]]$X)
prior_parameters <- list(
  a = 2,
  b = 3,
  tilde_a = 1/2,
  tilde_b = 1/2,
  alpha = rep(1/(K+1), K+1),
  mu0 = rep(0, p),
  V0 = diag(p)
)

# Get prior parameters
a <- prior_parameters$a
b <- prior_parameters$b
tilde_a <- prior_parameters$tilde_a
tilde_b <- prior_parameters$tilde_b
alpha <- prior_parameters$alpha
mu0 <- prior_parameters$mu0
V0 <- prior_parameters$V0

data <- sce1_data[[1]]

# fit lm
data_lm <- data.frame(y = data$y)
data_lm <- cbind(data_lm, data$X[, 2:ncol(data$X)])  # Exclude intercept

fit <- lm(y ~ ., data = data_lm)
summary(fit)
beta_mle <- coef(fit)

# Prepare data for Stan
stan_data_npp <- list(
  K = nrow(data$dims_X0),
  n0 =data$dims_X0[,1],
  n = nrow(data$X),
  p = ncol(data$X),
  len_X0 = length(data$X0_flat),
  len_y0 = length(data$y0_flat),
  X0 = data$X0_flat,
  startid_X0 = data$startid_X0,
  y0 = data$y0_flat,
  startid_y0 = data$startid_y0,
  X = data$X,
  y = as.vector(data$y),
  a = a,
  b = b,
  V0 = V0,
  mu0 = mu0,
  tilde_a = tilde_a,
  tilde_b = tilde_b,
  post = 1,
  seq = 0
)
stan_data_nppseq <- list(
  K = nrow(data$dims_X0),
  n0 =data$dims_X0[,1],
  n = nrow(data$X),
  p = ncol(data$X),
  len_X0 = length(data$X0_flat),
  len_y0 = length(data$y0_flat),
  X0 = data$X0_flat,
  startid_X0 = data$startid_X0,
  y0 = data$y0_flat,
  startid_y0 = data$startid_y0,
  X = data$X,
  y = as.vector(data$y),
  a = a,
  b = b,
  V0 = V0,
  mu0 = mu0,
  tilde_a = tilde_a,
  tilde_b = tilde_b,
  post = 1,
  seq = 1
)
stan_data_onpp <- list(
  K = nrow(data$dims_X0),
  n0 =data$dims_X0[,1],
  n = nrow(data$X),
  p = ncol(data$X),
  len_X0 = length(data$X0_flat),
  len_y0 = length(data$y0_flat),
  X0 = data$X0_flat,
  startid_X0 = data$startid_X0,
  y0 = data$y0_flat,
  startid_y0 = data$startid_y0,
  X = data$X,
  y = as.vector(data$y),
  a = a,
  b = b,
  V0 = V0,
  mu0 = mu0,
  alpha = alpha,
  post = 1,
  seq = 0
)
stan_data_onppseq <- list(
  K = nrow(data$dims_X0),
  n0 =data$dims_X0[,1],
  n = nrow(data$X),
  p = ncol(data$X),
  len_X0 = length(data$X0_flat),
  len_y0 = length(data$y0_flat),
  X0 = data$X0_flat,
  startid_X0 = data$startid_X0,
  y0 = data$y0_flat,
  startid_y0 = data$startid_y0,
  X = data$X,
  y = as.vector(data$y),
  a = a,
  b = b,
  V0 = V0,
  mu0 = mu0,
  alpha = alpha,
  post = 1,
  seq = 1
)

# Compile the model
gamma_model <- cmdstan_model("code/models/gamma_lm.stan")
delta_model <- cmdstan_model("code/models/delta_lm.stan")

sample_delta_npp <- delta_model$sample(
  data = stan_data_npp,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 4000,
  iter_sampling = 2000,
  refresh = 0,
  adapt_delta = 0.99,
  step_size = 0.03
)
sample_delta_nppseq <- delta_model$sample(
  data = stan_data_nppseq,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 4000,
  iter_sampling = 2000,
  refresh = 0,
  adapt_delta = 0.99,
  step_size = 0.03
)
sample_delta_onpp <- gamma_model$sample(
  data = stan_data_onpp,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 4000,
  iter_sampling = 2000,
  refresh = 0,
  adapt_delta = 0.99,
  step_size = 0.03
)

sample_delta_onppseq <- gamma_model$sample(
  data = stan_data_onppseq,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 4000,
  iter_sampling = 2000,
  refresh = 0,
  adapt_delta = 0.99,
  step_size = 0.03
)

sample_delta_npp <- stan(file = "code/models/delta_lm.stan", 
                         data = stan_data_npp,
                         chains = 4,
                         cores = 4,
                         warmup = 4000,
                         iter = 6000)

summary(sample_delta_npp)$summary

sample_delta_nppseq <- stan(file = "code/models/delta_lm.stan", 
                            data = stan_data_nppseq,
                            chains = 4,
                            cores = 4,
                            warmup = 4000,
                            iter = 6000)

summary(sample_delta_nppseq)$summary

sample_delta_onpp <- stan(file = "code/models/gamma_lm.stan", 
                         data = stan_data_onpp,
                         chains = 4,
                         cores = 4,
                         warmup = 4000,
                         iter = 6000)

summary(sample_delta_onpp)$summary

sample_delta_onppseq <- stan(file = "code/models/gamma_lm.stan", 
                            data = stan_data_onppseq,
                            chains = 4,
                            cores = 4,
                            warmup = 4000,
                            iter = 6000)

summary(sample_delta_onppseq)$summary


draws_delta_npp <- extract(sample_delta_npp)$delta %>%
  as_draws_matrix()
draws_delta_nppseq <- extract(sample_delta_nppseq)$delta %>%
  as_draws_matrix()
draws_delta_onpp <- extract(sample_delta_onpp)$delta %>%
  as_draws_matrix()
draws_delta_onppseq <- extract(sample_delta_onppseq)$delta %>%
  as_draws_matrix()
draws_beta_npp <- extract(sample_delta_npp)$beta %>%
  as_draws_matrix()
draws_beta_nppseq <- extract(sample_delta_nppseq)$beta %>%
  as_draws_matrix()
draws_beta_onpp <- extract(sample_delta_onpp)$beta %>%
  as_draws_matrix()
draws_beta_onppseq <- extract(sample_delta_onppseq)$beta %>%
  as_draws_matrix()
draws_sigma_npp <- matrix(extract(sample_delta_npp)$sigma, ncol = 1)
draws_sigma_nppseq <- matrix(extract(sample_delta_nppseq)$sigma, ncol = 1)
draws_sigma_onpp <- matrix(extract(sample_delta_onpp)$sigma, ncol = 1)
draws_sigma_onppseq <- matrix(extract(sample_delta_onppseq)$sigma, ncol = 1)

delta1 <- cbind(draws_delta_npp[,1], 
                draws_delta_nppseq[,1],
                draws_delta_onpp[,1],
                draws_delta_onppseq[,1])
delta2 <- cbind(draws_delta_npp[,2],
                draws_delta_nppseq[,2],
                draws_delta_onpp[,2],
                draws_delta_onppseq[,2])
delta3 <- cbind(draws_delta_npp[,3],
                draws_delta_nppseq[,3],
                draws_delta_onpp[,3],
                draws_delta_onppseq[,3])

beta1 <- cbind(draws_beta_npp[,1],
               draws_beta_nppseq[,1],
               draws_beta_onpp[,1],
               draws_beta_onppseq[,1])
beta2 <- cbind(draws_beta_npp[,2],
               draws_beta_nppseq[,2],
               draws_beta_onpp[,2],
               draws_beta_onppseq[,2])
beta3 <- cbind(draws_beta_npp[,3],
               draws_beta_nppseq[,3],
               draws_beta_onpp[,3],
               draws_beta_onppseq[,3])
sigma <- cbind(draws_sigma_npp[,1],
               draws_sigma_nppseq[,1],
               draws_sigma_onpp[,1],
               draws_sigma_onppseq[,1])

# Plotting
box_delta1 <- plot_boxplot(delta1, expression(delta[1]))
box_delta2 <- plot_boxplot(delta2, expression(delta[2]))
box_delta3 <- plot_boxplot(delta3, expression(delta[3]))
box_beta1 <- plot_boxplot(beta1, expression(beta[1])) + 
  geom_hline(yintercept = data$beta[1],
             color = "red", linetype = "dashed", size = 1) +
  geom_hline(yintercept = beta_mle[1],
             color = "blue", linetype = "dashed", size = 1)
box_beta2 <- plot_boxplot(beta2, expression(beta[2])) +
  geom_hline(yintercept = data$beta[2],
             color = "red", linetype = "dashed", size = 1) +
  geom_hline(yintercept = beta_mle[2],
             color = "blue", linetype = "dashed", size = 1)
box_beta3 <- plot_boxplot(beta3, expression(beta[3])) +
  geom_hline(yintercept = data$beta[3],
             color = "red", linetype = "dashed", size = 1) +
  geom_hline(yintercept = beta_mle[3],
             color = "blue", linetype = "dashed", size = 1)
box_sigma <- plot_boxplot(sigma, expression(sigma^2)) +
  geom_hline(yintercept = data$sg,
             color = "red", linetype = "dashed", size = 1)

box_delta1
box_delta2
box_delta3
box_beta1
box_beta2
box_beta3
box_sigma

# Shinystan diagnostics
my_sso <- launch_shinystan(sample_delta_onpp)

pairs(sample_delta_onpp, pars = c("delta", "beta", "sigma"), las = 1)

traceplot(sample_delta_onpp, pars = c("gamma","delta"), nrow = 2)

#-------------------------------------------------------------------------------

y <- data.frame(m=1/rgamma(1000, shape = 2, rate = 2))

ggplot(y, aes(x = m)) +
  geom_histogram(aes(y = ..density..), bins = 60, fill = "steelblue", alpha = 0.5) +
  geom_density(color = "red", size = 1) +
  labs(title = "Gamma Distribution",
       x = "x",
       y = "Density") +
  theme_minimal()+
  xlim(0, 2)

