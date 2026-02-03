# Load data
sim1.1 <- readRDS("results/sim_data/lm/sce_I_I.rds")
sim1.2 <- readRDS("results/sim_data/lm/sce_I_II.rds")
sim1.3 <- readRDS("results/sim_data/lm/sce_I_III.rds")
sim2.1 <- readRDS("results/sim_data/lm/sce_II_I.rds")
sim2.2 <- readRDS("results/sim_data/lm/sce_II_II.rds")
sim2.3 <- readRDS("results/sim_data/lm/sce_II_III.rds")
sim3.1 <- readRDS("results/sim_data/lm/sce_III_I.rds")
sim3.2 <- readRDS("results/sim_data/lm/sce_III_II.rds")
sim3.3 <- readRDS("results/sim_data/lm/sce_III_III.rds")
sim4.1 <- readRDS("results/sim_data/lm/sce_IV_I.rds")
sim4.2 <- readRDS("results/sim_data/lm/sce_IV_II.rds")
sim4.3 <- readRDS("results/sim_data/lm/sce_IV_III.rds")

sim_data_list <- list(
  sim1.1, sim1.2, sim1.3,
  sim2.1, sim2.2, sim2.3,
  sim3.1, sim3.2, sim3.3,
  sim4.1, sim4.2, sim4.3
)

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
    # parallel_chains = 4,
    iter_warmup = 4000,
    iter_sampling = 2000,
    refresh = 0
    # adapt_delta = 0.99,
    # step_size = 0.03
  )
  sample_delta_nppseq <- delta_model$sample(
    data = stan_data_nppseq,
    chains = 4,
    # parallel_chains = 4,
    iter_warmup = 4000,
    iter_sampling = 2000,
    refresh = 0
    # adapt_delta = 0.99,
    # step_size = 0.03
  )
  sample_delta_onpp <- gamma_model$sample(
    data = stan_data_onpp,
    chains = 4,
    # parallel_chains = 4,
    iter_warmup = 4000,
    iter_sampling = 2000,
    refresh = 0,
    adapt_delta = 0.99
    # step_size = 0.03
  )
  
  sample_delta_onppseq <- gamma_model$sample(
    data = stan_data_onppseq,
    chains = 4,
    # parallel_chains = 4,
    iter_warmup = 4000,
    iter_sampling = 2000,
    refresh = 0,
    adapt_delta = 0.99
    # step_size = 0.03
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
  
  
  divergences_npp <- sum(((sample_delta_npp$sampler_diagnostics())[, , "divergent__"] %>% 
                            as_draws_df())$divergent__)
  divergences_nppseq <- sum(((sample_delta_nppseq$sampler_diagnostics())[, , "divergent__"] %>% 
                               as_draws_df())$divergent__)
  divergences_onpp <- sum(((sample_delta_onpp$sampler_diagnostics())[, , "divergent__"] %>%
                             as_draws_df())$divergent__)
  divergences_onppseq <- sum(((sample_delta_onppseq$sampler_diagnostics())[, , "divergent__"] %>%
                                as_draws_df())$divergent__)
  
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
    draws_sigma_onppseq = draws_sigma_onppseq,
    divergences = c(divergences_npp,
                     divergences_nppseq,
                     divergences_onpp,
                     divergences_onppseq)
  ))
}


# define prior parameters
K <- length(sim1.1[[1]]$X0)
p <- ncol(sim1.1[[1]]$X)
hyper <- list(
  a = 2,
  b = 1,
  tilde_a = 1/2,
  tilde_b = 1/2,
  alpha = rep(1/(K+1), K+1),
  mu0 = rep(0, p),
  V0 = diag(p)
)

# Fitting datasets in parallel
samples_sces <- lapply(seq_along(sim_data_list), function(i) {
  sce_data <- sim_data_list[[i]]
  # samples <- mclapply(sce_data, function(data) {
  #   fit_lm(data, 
  #          list(gamma_model = gamma_model,
  #               delta_model = delta_model),
  #          hyper)
  # },
  # mc.cores = 10)
  # create cluster with N workers
  cl <- makeCluster(15, type = "PSOCK")
  
  # export needed objects/functions to each worker
  clusterExport(cl, c("fit_lm", "gamma_model", "delta_model", "hyper"))
  clusterEvalQ(cl, {
    library(cmdstanr)
    library(posterior)
    library(dplyr)
  })
  
  samples <- parLapply(
    cl,
    sce_data,
    function(data) {
      fit_lm(
        data,
        list(
          gamma_model = gamma_model,
          delta_model = delta_model
        ),
        hyper
      )
    }
  )
  
  stopCluster(cl)
  sce <- as.roman(ceiling(i/3))
  sce <- ifelse(i %% 3 == 1, paste0(sce,"_I"), ifelse(i %% 3 == 2, paste0(sce,"_II"), paste0(sce,"_III")))
  saveRDS(samples, file = paste0("results/samples/lm/sce_", sce, ".rds"))
})

