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
delta_model_bin_os <- cmdstan_model("code/models/delta_bin_os.stan")

#-------------------------------------------------------------------------------
# Run simulations for the binomial model with different scenarios
#-------------------------------------------------------------------------------

# define scenarios
## compatible order
### compatible order and high congruence
sce1.1 <- list(n0 = c(30, 30, 30),
              n = 90,
              a = 1/2,
              b = 1/2,
              al = 1,
              bl = 1,
              theta0 = c(0.3, 0.4, 0.5),
              theta = 0.5,
              alpha = rep(1/4, 4)
)
### compatible order and small congruence
sce1.2 <- list(n0 = c(30, 30, 30),
             n = 90,
             a = 1/2,
             b = 1/2,
             al = 1,
             bl = 1,
             theta0 = c(0.2, 0.3, 0.4),
             theta = 0.5,
             alpha = rep(1/4, 4)
)
### compatible order and no congruence 
sce1.3 <- list(n0 = c(30, 30, 30),
             n = 90,
             a = 1/2,
             b = 1/2,
             al = 1,
             bl = 1,
             theta0 = c(0.05, 0.15, 0.25),
             theta = 0.5,
             alpha = rep(1/4, 4)
)

## almost compatible order
### almost compatible order and high congruence
sce2.1 <- list(n0 = c(30, 30, 30),
             n = 90,
             a = 1/2,
             b = 1/2,
             al = 1,
             bl = 1,
             theta0 = c(0.3, 0.5, 0.4),
             theta = 0.5,
             alpha = rep(1/4, 4)
)
### almost compatible order and small congruence
sce2.2 <- list(n0 = c(30, 30, 30),
             n = 90,
             a = 1/2,
             b = 1/2,
             al = 1,
             bl = 1,
             theta0 = c(0.2, 0.4, 0.3),
             theta = 0.5,
             alpha = rep(1/4, 4)
)
### almost compatible order and no congruence
sce2.3 <- list(n0 = c(30, 30, 30),
             n = 90,
             a = 1/2,
             b = 1/2,
             al = 1,
             bl = 1,
             theta0 = c(0.05, 0.25, 0.15),
             theta = 0.5,
             alpha = rep(1/4, 4)
)

## incompatible order
### incompatible order and huge congruence
sce3.1 <- list(n0 = c(30, 30, 30),
             n = 90,
             a = 1/2,
             b = 1/2,
             al = 1,
             bl = 1,
             theta0 = c(0.5, 0.4, 0.3),
             theta = 0.5,
             alpha = rep(1/4, 4)
)
### incompatible order and small congruence
sce3.2 <- list(n0 = c(30, 30, 30),
             n = 90,
             a = 1/2,
             b = 1/2,
             al = 1,
             bl = 1,
             theta0 = c(0.4, 0.3, 0.2),
             theta = 0.5,
             alpha = rep(1/4, 4)
)
### incompatible order and no congruence
sce3.3 <- list(n0 = c(30, 30, 30),
             n = 90,
             a = 1/2,
             b = 1/2,
             al = 1,
             bl = 1,
             theta0 = c(0.25, 0.15, 0.05),
             theta = 0.5,
             alpha = rep(1/4, 4)
)

## neutral order
### neutral order and huge congruence
sce4.1 <- list(n0 = c(30, 30, 30),
             n = 90,
             a = 1/2,
             b = 1/2,
             al = 1,
             bl = 1,
             theta0 = c(0.5, 0.5, 0.5),
             theta = 0.5,
             alpha = rep(1/4, 4)
)
### neutral order and small congruence
sce4.2 <- list(n0 = c(30, 30, 30),
             n = 90,
             a = 1/2,
             b = 1/2,
             al = 1,
             bl = 1,
             theta0 = c(0.3, 0.3, 0.3),
             theta = 0.5,
             alpha = rep(1/4, 4)
)
### neutral order and no congruence
sce4.3 <- list(n0 = c(30, 30, 30),
             n = 90,
             a = 1/2,
             b = 1/2,
             al = 1,
             bl = 1,
             theta0 = c(0.1, 0.1, 0.1),
             theta = 0.5,
             alpha = rep(1/4, 4)
)

# run simulations
sim1.1 <- sim_sce(model = "bin",14, 200, sce1.1, gamma_model_bin, delta_model_bin)
sim1.2 <- sim_sce(model = "bin",14, 200, sce1.2, gamma_model_bin, delta_model_bin)
sim1.3 <- sim_sce(model = "bin",14, 200, sce1.3, gamma_model_bin, delta_model_bin)
sim2.1 <- sim_sce(model = "bin",14, 200, sce2.1, gamma_model_bin, delta_model_bin)
sim2.2 <- sim_sce(model = "bin",14, 200, sce2.2, gamma_model_bin, delta_model_bin)
sim2.3 <- sim_sce(model = "bin",14, 200, sce2.3, gamma_model_bin, delta_model_bin)
sim3.1 <- sim_sce(model = "bin",14, 200, sce3.1, gamma_model_bin, delta_model_bin)
sim3.2 <- sim_sce(model = "bin",14, 200, sce3.2, gamma_model_bin, delta_model_bin)
sim3.3 <- sim_sce(model = "bin",14, 200, sce3.3, gamma_model_bin, delta_model_bin)
sim4.1 <- sim_sce(model = "bin",14, 200, sce4.1, gamma_model_bin, delta_model_bin)
sim4.2 <- sim_sce(model = "bin",14, 200, sce4.2, gamma_model_bin, delta_model_bin)
sim4.3 <- sim_sce(model = "bin",14, 200, sce4.3, gamma_model_bin, delta_model_bin)



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

#-------------------------------------------------------------------------------
# Run simulations for the binomial model with order strongly suggested
#-------------------------------------------------------------------------------

### incompatible order and huge congruence
sce3.1_os <- list(n0 = c(30, 30, 30),
               n = 30,
               a = 1/2,
               b = 1/2,
               al = c(1, 10, 10),
               bl = c(10, 10, 1),
               theta0 = c(0.5, 0.4, 0.3),
               theta = 0.5,
               alpha = rep(1/4, 4)
)
### incompatible order and small congruence
sce3.2_os <- list(n0 = c(30, 30, 30),
               n = 30,
               a = 1/2,
               b = 1/2,
               al = c(1, 10, 10),
               bl = c(10, 10, 1),
               theta0 = c(0.4, 0.3, 0.2),
               theta = 0.5,
               alpha = rep(1/4, 4)
)

### incompatible order and no congruence
sce3.3_os <- list(n0 = c(30, 30, 30),
               n = 30,
               a = 1/2,
               b = 1/2,
               al = c(1, 10, 10),
               bl = c(10, 10, 1),
               theta0 = c(0.25, 0.15, 0.05),
               theta = 0.5,
               alpha = rep(1/4, 4)
)

sim3.1_os <- sim_sce(model = "bin",14, 200, sce3.1_os, gamma_model_bin, delta_model_bin_os)
sim3.2_os <- sim_sce(model = "bin",14, 200, sce3.2_os, gamma_model_bin, delta_model_bin_os)
sim3.3_os <- sim_sce(model = "bin",14, 200, sce3.3_os, gamma_model_bin, delta_model_bin_os)

save(sim3.1_os, file = "results/samples/bin/sim3_1_os_bin.RData")
save(sim3.2_os, file = "results/samples/bin/sim3_2_os_bin.RData")
save(sim3.3_os, file = "results/samples/bin/sim3_3_os_bin.RData")
