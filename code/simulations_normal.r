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
gamma_model_normal <- cmdstan_model("code/models/normal_fixed_var_onpp.stan")
delta_model_normal <- cmdstan_model("code/models/normal_fixed_var_npp.stan")

#-------------------------------------------------------------------------------
# Run simulations for the normal model with fixed variance
#-------------------------------------------------------------------------------

# define scenarios
## compatible order
### compatible order and high congruence
sce1.1 <- list(n0 = c(100, 100, 100),
              n = 100,
              sigma0 = 1,
              sigmah = c(1, 1, 1),
              sigma = 1,
              mu0 = 0,
              al = 1/2,
              bl = 1/2,
              theta0 = c(4.6, 4.8, 5),
              theta = 5,
              alpha = rep(1/4, 4)
)
### compatible order and small congruence
sce1.2 <- list(n0 = c(100, 100, 100),
              n = 100,
              sigma0 = 1,
              sigmah = c(1, 1, 1),
              sigma = 1,
              mu0 = 0,
              al = 1/2,
              bl = 1/2,
              theta0 = c(3.6, 3.8, 4),
              theta = 5,
              alpha = rep(1/4, 4)
)
### compatible order and no congruence 
sce1.3 <- list(n0 = c(100, 100, 100),
              n = 100,
              sigma0 = 1,
              sigmah = c(1, 1, 1),
              sigma = 1,
              mu0 = 0,
              al = 1/2,
              bl = 1/2,
              theta0 = c(2.6, 2.8, 3),
              theta = 5,
              alpha = rep(1/4, 4)
)

## almost compatible order
### almost compatible order and high congruence
sce2.1 <- list(n0 = c(100, 100, 100),
              n = 100,
              sigma0 = 1,
              sigmah = c(1, 1, 1),
              sigma = 1,
              mu0 = 0,
              al = 1/2,
              bl = 1/2,
              theta0 = c(4.6, 5, 4.8),
              theta = 5,
              alpha = rep(1/4, 4)
)
### almost compatible order and small congruence
sce2.2 <- list(n0 = c(100, 100, 100),
              n = 100,
              sigma0 = 1,
              sigmah = c(1, 1, 1),
              sigma = 1,
              mu0 = 0,
              al = 1/2,
              bl = 1/2,
              theta0 = c(3.6, 4, 3.8),
              theta = 5,
              alpha = rep(1/4, 4)
)
### almost compatible order and no congruence
sce2.3 <- list(n0 = c(100, 100, 100),
              n = 100,
              sigma0 = 1,
              sigmah = c(1, 1, 1),
              sigma = 1,
              mu0 = 0,
              al = 1/2,
              bl = 1/2,
              theta0 = c(2.6, 3, 2.8),
              theta = 5,
              alpha = rep(1/4, 4)
)

## incompatible order
### incompatible order and huge congruence
sce3.1 <- list(n0 = c(100, 100, 100),
              n = 100,
              sigma0 = 1,
              sigmah = c(1, 1, 1),
              sigma = 1,
              mu0 = 0,
              al = 1/2,
              bl = 1/2,
              theta0 = c(5, 4.8, 4.6),
              theta = 5,
              alpha = rep(1/4, 4)
)
### incompatible order and small congruence
sce3.2 <- list(n0 = c(100, 100, 100),
              n = 100,
              sigma0 = 1,
              sigmah = c(1, 1, 1),
              sigma = 1,
              mu0 = 0,
              al = 1/2,
              bl = 1/2,
              theta0 = c(4, 3.8, 3.6),
              theta = 5,
              alpha = rep(1/4, 4)
)
### incompatible order and no congruence
sce3.3 <- list(n0 = c(100, 100, 100),
              n = 100,
              sigma0 = 1,
              sigmah = c(1, 1, 1),
              sigma = 1,
              mu0 = 0,
              al = 1/2,
              bl = 1/2,
              theta0 = c(3, 2.8, 2.6),
              theta = 5,
              alpha = rep(1/4, 4)
)

## neutral order
### neutral order and huge congruence
sce4.1 <- list(n0 = c(100, 100, 100),
              n = 100,
              sigma0 = 1,
              sigmah = c(1, 1, 1),
              sigma = 1,
              mu0 = 0,
              al = 1/2,
              bl = 1/2,
              theta0 = c(5, 5, 5),
              theta = 5,
              alpha = rep(1/4, 4)
)
### neutral order and small congruence
sce4.2 <- list(n0 = c(100, 100, 100),
              n = 100,
              sigma0 = 1,
              sigmah = c(1, 1, 1),
              sigma = 1,
              mu0 = 0,
              al = 1/2,
              bl = 1/2,
              theta0 = c(4, 4, 4),
              theta = 5,
              alpha = rep(1/4, 4)
)
### neutral order and no congruence
sce4.3 <- list(n0 = c(100, 100, 100),
              n = 100,
              sigma0 = 1,
              sigmah = c(1, 1, 1),
              sigma = 1,
              mu0 = 0,
              al = 1/2,
              bl = 1/2,
              theta0 = c(3, 3, 3),
              theta = 5,
              alpha = rep(1/4, 4)
)
## cancellation effect
sce5.1 <- list(n0 = c(100, 100, 100),
              n = 100,
              sigma0 = 1,
              sigmah = c(1, 1, 1),
              sigma = 1,
              mu0 = 0,
              al = 1/2,
              bl = 1/2,
              theta0 = c(3, 7, 5),
              theta = 5,
              alpha = rep(1/4, 4)
)

# run simulations
n_cores <- 15
n_sim <- 200
sim1.1 <- sim_sce(model = "normal_fixed_var", n_cores,  n_sim, sce1.1, gamma_model_normal, delta_model_normal)
sim1.2 <- sim_sce(model = "normal_fixed_var", n_cores,  n_sim, sce1.2, gamma_model_normal, delta_model_normal)
sim1.3 <- sim_sce(model = "normal_fixed_var", n_cores,  n_sim, sce1.3, gamma_model_normal, delta_model_normal)
sim2.1 <- sim_sce(model = "normal_fixed_var", n_cores,  n_sim, sce2.1, gamma_model_normal, delta_model_normal)
sim2.2 <- sim_sce(model = "normal_fixed_var", n_cores,  n_sim, sce2.2, gamma_model_normal, delta_model_normal)
sim2.3 <- sim_sce(model = "normal_fixed_var", n_cores,  n_sim, sce2.3, gamma_model_normal, delta_model_normal)
sim3.1 <- sim_sce(model = "normal_fixed_var", n_cores,  n_sim, sce3.1, gamma_model_normal, delta_model_normal)
sim3.2 <- sim_sce(model = "normal_fixed_var", n_cores,  n_sim, sce3.2, gamma_model_normal, delta_model_normal)
sim3.3 <- sim_sce(model = "normal_fixed_var", n_cores,  n_sim, sce3.3, gamma_model_normal, delta_model_normal)
sim4.1 <- sim_sce(model = "normal_fixed_var", n_cores,  n_sim, sce4.1, gamma_model_normal, delta_model_normal)
sim4.2 <- sim_sce(model = "normal_fixed_var", n_cores,  n_sim, sce4.2, gamma_model_normal, delta_model_normal)
sim4.3 <- sim_sce(model = "normal_fixed_var", n_cores,  n_sim, sce4.3, gamma_model_normal, delta_model_normal)
sim5.1 <- sim_sce(model = "normal_fixed_var", n_cores,  n_sim, sce5.1, gamma_model_normal, delta_model_normal)

# save results
qs2::qs_save(sim1.1, file = "results/samples/normal/sim1_1_normal_fixed_var.qs2")
qs2::qs_save(sim1.2, file = "results/samples/normal/sim1_2_normal_fixed_var.qs2")
qs2::qs_save(sim1.3, file = "results/samples/normal/sim1_3_normal_fixed_var.qs2")
qs2::qs_save(sim2.1, file = "results/samples/normal/sim2_1_normal_fixed_var.qs2")
qs2::qs_save(sim2.2, file = "results/samples/normal/sim2_2_normal_fixed_var.qs2")
qs2::qs_save(sim2.3, file = "results/samples/normal/sim2_3_normal_fixed_var.qs2")
qs2::qs_save(sim3.1, file = "results/samples/normal/sim3_1_normal_fixed_var.qs2")
qs2::qs_save(sim3.2, file = "results/samples/normal/sim3_2_normal_fixed_var.qs2")
qs2::qs_save(sim3.3, file = "results/samples/normal/sim3_3_normal_fixed_var.qs2")
qs2::qs_save(sim4.1, file = "results/samples/normal/sim4_1_normal_fixed_var.qs2")
qs2::qs_save(sim4.2, file = "results/samples/normal/sim4_2_normal_fixed_var.qs2")
qs2::qs_save(sim4.3, file = "results/samples/normal/sim4_3_normal_fixed_var.qs2")
qs2::qs_save(sim5.1, file = "results/samples/normal/sim5_1_normal_fixed_var.qs2")
