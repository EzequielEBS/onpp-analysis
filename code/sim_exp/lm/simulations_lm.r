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
library(qs2)

# run auxiliary functions
source("code/aux_fun_sim.R")

gamma_model_lm <- cmdstan_model("code/models/lm/gamma_lm.stan")
delta_model_lm <- cmdstan_model("code/models/lm/delta_lm.stan")

sce1.1 <- qs_read("results/sim_data/lm/sceI_I.qs2")
sce1.2 <- qs_read("results/sim_data/lm/sceI_II.qs2")
sce1.3 <- qs_read("results/sim_data/lm/sceI_III.qs2")
sce2.1 <- qs_read("results/sim_data/lm/sceII_I.qs2")
sce2.2 <- qs_read("results/sim_data/lm/sceII_II.qs2")
sce2.3 <- qs_read("results/sim_data/lm/sceII_III.qs2")
sce3.1 <- qs_read("results/sim_data/lm/sceIII_I.qs2")
sce3.2 <- qs_read("results/sim_data/lm/sceIII_II.qs2")
sce3.3 <- qs_read("results/sim_data/lm/sceIII_III.qs2")
sce4.1 <- qs_read("results/sim_data/lm/sceIV_I.qs2")
sce4.2 <- qs_read("results/sim_data/lm/sceIV_II.qs2")
sce4.3 <- qs_read("results/sim_data/lm/sceIV_III.qs2")

# run simulations
n_cores <- 15
n_sim <- 200
sim1.1 <- sim_sce(model = "lm", n_cores,  n_sim, sce1.1, gamma_model_lm, delta_model_lm)
qs_save(sim1.1,
  file = "results/samples/lm/sim1_1_lm.qs2"
)

sim1.2 <- sim_sce(model = "lm", n_cores,  n_sim, sce1.2, gamma_model_lm, delta_model_lm)
qs_save(sim1.2,
  file = "results/samples/lm/sim1_2_lm.qs2"
)

sim1.3 <- sim_sce(model = "lm", n_cores,  n_sim, sce1.3, gamma_model_lm, delta_model_lm)
qs_save(sim1.3,
  file = "results/samples/lm/sim1_3_lm.qs2"
)

sim2.1 <- sim_sce(model = "lm", n_cores,  n_sim, sce2.1, gamma_model_lm, delta_model_lm)
qs_save(sim2.1,
  file = "results/samples/lm/sim2_1_lm.qs2"
)

sim2.2 <- sim_sce(model = "lm", n_cores,  n_sim, sce2.2, gamma_model_lm, delta_model_lm)
qs_save(sim2.2,
  file = "results/samples/lm/sim2_2_lm.qs2"
)

sim2.3 <- sim_sce(model = "lm", n_cores,  n_sim, sce2.3, gamma_model_lm, delta_model_lm)
qs_save(sim2.3,
  file = "results/samples/lm/sim2_3_lm.qs2"
)

sim3.1 <- sim_sce(model = "lm", n_cores,  n_sim, sce3.1, gamma_model_lm, delta_model_lm)
qs_save(sim3.1,
  file = "results/samples/lm/sim3_1_lm.qs2"
)

sim3.2 <- sim_sce(model = "lm", n_cores,  n_sim, sce3.2, gamma_model_lm, delta_model_lm)
qs_save(sim3.2,
  file = "results/samples/lm/sim3_2_lm.qs2"
)

sim3.3 <- sim_sce(model = "lm", n_cores,  n_sim, sce3.3, gamma_model_lm, delta_model_lm)
qs_save(sim3.3,
  file = "results/samples/lm/sim3_3_lm.qs2"
)

sim4.1 <- sim_sce(model = "lm", n_cores,  n_sim, sce4.1, gamma_model_lm, delta_model_lm)
qs_save(sim4.1,
  file = "results/samples/lm/sim4_1_lm.qs2"
)

sim4.2 <- sim_sce(model = "lm", n_cores,  n_sim, sce4.2, gamma_model_lm, delta_model_lm)
qs_save(sim4.2,
  file = "results/samples/lm/sim4_2_lm.qs2"
)

sim4.3 <- sim_sce(model = "lm", n_cores,  n_sim, sce4.3, gamma_model_lm, delta_model_lm)
qs_save(sim4.3,
  file = "results/samples/lm/sim4_3_lm.qs2"
)
