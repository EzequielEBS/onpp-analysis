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

# run auxiliary functions
source("code/aux_fun_sim.R")

# Compile the model
gamma_model_bin <- cmdstan_model("code/models/gamma_bin.stan")
delta_model_bin <- cmdstan_model("code/models/delta_bin.stan")

#-------------------------------------------------------------------------------
# Run simulations for the binomial model with different scenarios
#-------------------------------------------------------------------------------

# define scenarios
## compatible order
### compatible order and high congruence
sce1.1 <- list(n0 = c(30, 30, 30),
              n = 30,
              a = 1/2,
              b = 1/2,
              al = 1/2,
              bl = 1/2,
              theta0 = c(0.3, 0.4, 0.5),
              theta = 0.5,
              alpha = rep(1/4, 4)
)
### compatible order and small congruence
sce1.2 <- list(n0 = c(30, 30, 30),
             n = 30,
             a = 1/2,
             b = 1/2,
             al = 1/2,
             bl = 1/2,
             theta0 = c(0.2, 0.3, 0.4),
             theta = 0.5,
             alpha = rep(1/4, 4)
)
### compatible order and no congruence 
sce1.3 <- list(n0 = c(30, 30, 30),
             n = 30,
             a = 1/2,
             b = 1/2,
             al = 1/2,
             bl = 1/2,
             theta0 = c(0.05, 0.15, 0.25),
             theta = 0.5,
             alpha = rep(1/4, 4)
)

## almost compatible order
### almost compatible order and high congruence
sce2.1 <- list(n0 = c(30, 30, 30),
             n = 30,
             a = 1/2,
             b = 1/2,
             al = 1/2,
             bl = 1/2,
             theta0 = c(0.3, 0.5, 0.4),
             theta = 0.5,
             alpha = rep(1/4, 4)
)
### almost compatible order and small congruence
sce2.2 <- list(n0 = c(30, 30, 30),
             n = 30,
             a = 1/2,
             b = 1/2,
             al = 1/2,
             bl = 1/2,
             theta0 = c(0.2, 0.4, 0.3),
             theta = 0.5,
             alpha = rep(1/4, 4)
)
### almost compatible order and no congruence
sce2.3 <- list(n0 = c(30, 30, 30),
             n = 30,
             a = 1/2,
             b = 1/2,
             al = 1/2,
             bl = 1/2,
             theta0 = c(0.05, 0.25, 0.15),
             theta = 0.5,
             alpha = rep(1/4, 4)
)

## incompatible order
### incompatible order and huge congruence
sce3.1 <- list(n0 = c(30, 30, 30),
             n = 30,
             a = 1/2,
             b = 1/2,
             al = 1/2,
             bl = 1/2,
             theta0 = c(0.5, 0.4, 0.3),
             theta = 0.5,
             alpha = rep(1/4, 4)
)
### incompatible order and small congruence
sce3.2 <- list(n0 = c(30, 30, 30),
             n = 30,
             a = 1/2,
             b = 1/2,
             al = 1/2,
             bl = 1/2,
             theta0 = c(0.4, 0.3, 0.2),
             theta = 0.5,
             alpha = rep(1/4, 4)
)
### incompatible order and no congruence
sce3.3 <- list(n0 = c(30, 30, 30),
             n = 30,
             a = 1/2,
             b = 1/2,
             al = 1/2,
             bl = 1/2,
             theta0 = c(0.25, 0.15, 0.05),
             theta = 0.5,
             alpha = rep(1/4, 4)
)

## neutral order
### neutral order and huge congruence
sce4.1 <- list(n0 = c(30, 30, 30),
             n = 30,
             a = 1/2,
             b = 1/2,
             al = 1/2,
             bl = 1/2,
             theta0 = c(0.5, 0.5, 0.5),
             theta = 0.5,
             alpha = rep(1/4, 4)
)
### neutral order and small congruence
sce4.2 <- list(n0 = c(30, 30, 30),
             n = 30,
             a = 1/2,
             b = 1/2,
             al = 1/2,
             bl = 1/2,
             theta0 = c(0.3, 0.3, 0.3),
             theta = 0.5,
             alpha = rep(1/4, 4)
)
### neutral order and no congruence
sce4.3 <- list(n0 = c(30, 30, 30),
             n = 30,
             a = 1/2,
             b = 1/2,
             al = 1/2,
             bl = 1/2,
             theta0 = c(0.1, 0.1, 0.1),
             theta = 0.5,
             alpha = rep(1/4, 4)
)
## cancelation effect
sce5.1 <- list(n0 = c(30, 30, 30),
             n = 30,
             a = 1/2,
             b = 1/2,
             al = 1/2,
             bl = 1/2,
             theta0 = c(0.1, 0.9, 0.5),
             theta = 0.5,
             alpha = rep(1/4, 4)
)

# run simulations
n_cores <- 15
n_sim <- 200
sim1.1 <- sim_sce(model = "bin", n_cores,  n_sim, sce1.1, gamma_model_bin, delta_model_bin)
sim1.2 <- sim_sce(model = "bin", n_cores,  n_sim, sce1.2, gamma_model_bin, delta_model_bin)
sim1.3 <- sim_sce(model = "bin", n_cores,  n_sim, sce1.3, gamma_model_bin, delta_model_bin)
sim2.1 <- sim_sce(model = "bin", n_cores,  n_sim, sce2.1, gamma_model_bin, delta_model_bin)
sim2.2 <- sim_sce(model = "bin", n_cores,  n_sim, sce2.2, gamma_model_bin, delta_model_bin)
sim2.3 <- sim_sce(model = "bin", n_cores,  n_sim, sce2.3, gamma_model_bin, delta_model_bin)
sim3.1 <- sim_sce(model = "bin", n_cores,  n_sim, sce3.1, gamma_model_bin, delta_model_bin)
sim3.2 <- sim_sce(model = "bin", n_cores,  n_sim, sce3.2, gamma_model_bin, delta_model_bin)
sim3.3 <- sim_sce(model = "bin", n_cores,  n_sim, sce3.3, gamma_model_bin, delta_model_bin)
sim4.1 <- sim_sce(model = "bin", n_cores,  n_sim, sce4.1, gamma_model_bin, delta_model_bin)
sim4.2 <- sim_sce(model = "bin", n_cores,  n_sim, sce4.2, gamma_model_bin, delta_model_bin)
sim4.3 <- sim_sce(model = "bin", n_cores,  n_sim, sce4.3, gamma_model_bin, delta_model_bin)
sim5.1 <- sim_sce(model = "bin", n_cores,  n_sim, sce5.1, gamma_model_bin, delta_model_bin)

samples_prior_sce3.1 <- sample_sce_bin(sce3.1, gamma_model_bin, delta_model_bin, post = 0)
samples_prior_sce3.2 <- sample_sce_bin(sce3.2, gamma_model_bin, delta_model_bin, post = 0)
samples_prior_sce3.3 <- sample_sce_bin(sce3.3, gamma_model_bin, delta_model_bin, post = 0)

# save results
save(sim1.1, file = "results/samples/bin/sim1_1_bin.RData")
save(sim1.2, file = "results/samples/bin/sim1_2_bin.RData")
save(sim1.3, file = "results/samples/bin/sim1_3_bin.RData")
save(sim2.1, file = "results/samples/bin/sim2_1_bin.RData")
save(sim2.2, file = "results/samples/bin/sim2_2_bin.RData")
save(sim2.3, file = "results/samples/bin/sim2_3_bin.RData")
save(sim3.1, file = "results/samples/bin/sim3_1_bin.RData")
save(sim3.2, file = "results/samples/bin/sim3_2_bin.RData")
save(sim3.3, file = "results/samples/bin/sim3_3_bin.RData")
save(sim4.1, file = "results/samples/bin/sim4_1_bin.RData")
save(sim4.2, file = "results/samples/bin/sim4_2_bin.RData")
save(sim4.3, file = "results/samples/bin/sim4_3_bin.RData")
save(sim5.1, file = "results/samples/bin/sim5_1_bin.RData")

save(samples_prior_sce3.1, file = "results/samples/bin/samples_prior_sce3_1_bin.RData")
save(samples_prior_sce3.2, file = "results/samples/bin/samples_prior_sce3_2_bin.RData")
save(samples_prior_sce3.3, file = "results/samples/bin/samples_prior_sce3_3_bin.RData")

# -------------------------------------------------------------------------------
# Flat prior for eta
## incompatible order
### incompatible order and huge congruence
sce3.1_eta_flat <- list(n0 = c(30, 30, 30),
             n = 30,
             a = 1/2,
             b = 1/2,
             al = 1,
             bl = 1,
             theta0 = c(0.5, 0.4, 0.3),
             theta = 0.5,
             alpha = rep(1/4, 4)
)
### incompatible order and small congruence
sce3.2_eta_flat <- list(n0 = c(30, 30, 30),
             n = 30,
             a = 1/2,
             b = 1/2,
             al = 1,
             bl = 1,
             theta0 = c(0.4, 0.3, 0.2),
             theta = 0.5,
             alpha = rep(1/4, 4)
)
### incompatible order and no congruence
sce3.3_eta_flat <- list(n0 = c(30, 30, 30),
             n = 30,
             a = 1/2,
             b = 1/2,
             al = 1,
             bl = 1,
             theta0 = c(0.25, 0.15, 0.05),
             theta = 0.5,
             alpha = rep(1/4, 4)
)

sim3.1_eta_flat <- sim_sce(model = "bin", n_cores,  n_sim, sce3.1_eta_flat, gamma_model_bin, delta_model_bin)
sim3.2_eta_flat <- sim_sce(model = "bin", n_cores,  n_sim, sce3.2_eta_flat, gamma_model_bin, delta_model_bin)
sim3.3_eta_flat <- sim_sce(model = "bin", n_cores,  n_sim, sce3.3_eta_flat, gamma_model_bin, delta_model_bin)

samples_prior_sce3.1_eta_flat <- sample_sce_bin(sce3.1_eta_flat, gamma_model_bin, delta_model_bin, post = 0)
samples_prior_sce3.2_eta_flat <- sample_sce_bin(sce3.2_eta_flat, gamma_model_bin, delta_model_bin, post = 0)
samples_prior_sce3.3_eta_flat <- sample_sce_bin(sce3.3_eta_flat, gamma_model_bin, delta_model_bin, post = 0)

save(sim3.1_eta_flat, file = "results/samples/bin/sim3_1_bin_eta_flat.RData")
save(sim3.2_eta_flat, file = "results/samples/bin/sim3_2_bin_eta_flat.RData")
save(sim3.3_eta_flat, file = "results/samples/bin/sim3_3_bin_eta_flat.RData")

save(samples_prior_sce3.1_eta_flat, file = "results/samples/bin/samples_prior_sce3_1_eta_flat_bin.RData")
save(samples_prior_sce3.2_eta_flat, file = "results/samples/bin/samples_prior_sce3_2_eta_flat_bin.RData")
save(samples_prior_sce3.3_eta_flat, file = "results/samples/bin/samples_prior_sce3_3_eta_flat_bin.RData")

# -------------------------------------------------------------------------------
# Flat prior for eta and gamma
## incompatible order
### incompatible order and huge congruence
sce3.1_eta_gamma_flat <- list(n0 = c(30, 30, 30),
             n = 30,
             a = 1/2,
             b = 1/2,
             al = 1,
             bl = 1,
             theta0 = c(0.5, 0.4, 0.3),
             theta = 0.5,
             alpha = rep(1, 4)
)

### incompatible order and small congruence
sce3.2_eta_gamma_flat <- list(n0 = c(30, 30, 30),
             n = 30,
             a = 1/2,
             b = 1/2,
             al = 1,
             bl = 1,
             theta0 = c(0.4, 0.3, 0.2),
             theta = 0.5,
             alpha = rep(1, 4)
)

### incompatible order and no congruence
sce3.3_eta_gamma_flat <- list(n0 = c(30, 30, 30),
             n = 30,
             a = 1/2,
             b = 1/2,
             al = 1,
             bl = 1,
             theta0 = c(0.25, 0.15, 0.05),
             theta = 0.5,
             alpha = rep(1, 4)
)

sim3.1_eta_gamma_flat <- sim_sce(model = "bin", n_cores,  n_sim, sce3.1_eta_gamma_flat, gamma_model_bin, delta_model_bin)
sim3.2_eta_gamma_flat <- sim_sce(model = "bin", n_cores,  n_sim, sce3.2_eta_gamma_flat, gamma_model_bin, delta_model_bin)
sim3.3_eta_gamma_flat <- sim_sce(model = "bin", n_cores,  n_sim, sce3.3_eta_gamma_flat, gamma_model_bin, delta_model_bin)

samples_prior_sce3.1_eta_gamma_flat <- sample_sce_bin(sce3.1_eta_gamma_flat, gamma_model_bin, delta_model_bin, post = 0)
samples_prior_sce3.2_eta_gamma_flat <- sample_sce_bin(sce3.2_eta_gamma_flat, gamma_model_bin, delta_model_bin, post = 0)
samples_prior_sce3.3_eta_gamma_flat <- sample_sce_bin(sce3.3_eta_gamma_flat, gamma_model_bin, delta_model_bin, post = 0)

save(sim3.1_eta_gamma_flat, file = "results/samples/bin/sim3_1_bin_eta_gamma_flat.RData")
save(sim3.2_eta_gamma_flat, file = "results/samples/bin/sim3_2_bin_eta_gamma_flat.RData")
save(sim3.3_eta_gamma_flat, file = "results/samples/bin/sim3_3_bin_eta_gamma_flat.RData")

save(samples_prior_sce3.1_eta_gamma_flat, file = "results/samples/bin/samples_prior_sce3_1_eta_gamma_flat_bin.RData")
save(samples_prior_sce3.2_eta_gamma_flat, file = "results/samples/bin/samples_prior_sce3_2_eta_gamma_flat_bin.RData")
save(samples_prior_sce3.3_eta_gamma_flat, file = "results/samples/bin/samples_prior_sce3_3_eta_gamma_flat_bin.RData")

# -------------------------------------------------------------------------------
# All flat priors
## incompatible order
### incompatible order and huge congruence
sce3.1_all_flat <- list(n0 = c(30, 30, 30),
             n = 30,
             a = 1,
             b = 1,
             al = 1,
             bl = 1,
             theta0 = c(0.5, 0.4, 0.3),
             theta = 0.5,
             alpha = rep(1, 4)
)

### incompatible order and small congruence
sce3.2_all_flat <- list(n0 = c(30, 30, 30),
             n = 30,
             a = 1,
             b = 1,
             al = 1,
             bl = 1,
             theta0 = c(0.4, 0.3, 0.2),
             theta = 0.5,
             alpha = rep(1, 4)
)

### incompatible order and no congruence
sce3.3_all_flat <- list(n0 = c(30, 30, 30),
             n = 30,
             a = 1,
             b = 1,
             al = 1,
             bl = 1,
             theta0 = c(0.25, 0.15, 0.05),
             theta = 0.5,
             alpha = rep(1, 4)
)

sim3.1_all_flat <- sim_sce(model = "bin", n_cores,  n_sim, sce3.1_all_flat, gamma_model_bin, delta_model_bin)
sim3.2_all_flat <- sim_sce(model = "bin", n_cores,  n_sim, sce3.2_all_flat, gamma_model_bin, delta_model_bin)
sim3.3_all_flat <- sim_sce(model = "bin", n_cores,  n_sim, sce3.3_all_flat, gamma_model_bin, delta_model_bin)

samples_prior_sce3.1_all_flat <- sample_sce_bin(sce3.1_all_flat, gamma_model_bin, delta_model_bin, post = 0)
samples_prior_sce3.2_all_flat <- sample_sce_bin(sce3.2_all_flat, gamma_model_bin, delta_model_bin, post = 0)
samples_prior_sce3.3_all_flat <- sample_sce_bin(sce3.3_all_flat, gamma_model_bin, delta_model_bin, post = 0)

save(sim3.1_all_flat, file = "results/samples/bin/sim3_1_bin_all_flat.RData")
save(sim3.2_all_flat, file = "results/samples/bin/sim3_2_bin_all_flat.RData")
save(sim3.3_all_flat, file = "results/samples/bin/sim3_3_bin_all_flat.RData")

save(samples_prior_sce3.1_all_flat, file = "results/samples/bin/samples_prior_sce3_1_all_flat_bin.RData")
save(samples_prior_sce3.2_all_flat, file = "results/samples/bin/samples_prior_sce3_2_all_flat_bin.RData")
save(samples_prior_sce3.3_all_flat, file = "results/samples/bin/samples_prior_sce3_3_all_flat_bin.RData")

