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
gamma_model_normal <- cmdstan_model("code/models/normal/normal_fixed_var_onpp.stan")
delta_model_normal <- cmdstan_model("code/models/normal/normal_fixed_var_npp.stan")

#-------------------------------------------------------------------------------
# Run simulations for the normal model with fixed variance
#-------------------------------------------------------------------------------
#
# All 14 scenarios (order I-IV x congruence I-III, plus the two
# "cancellation effect" scenarios V.1/V.2 that don't fit the order x
# congruence grid) share the exact same n0/n/sigma2/mu0/tau02/al/bl/theta/
# alpha and differ only in theta0. This used to be ~240 lines of
# copy-pasted `sceX.Y <- list(...)` / `simX.Y <- sim_sce(...)` /
# `qs_save(...)` blocks -- the same repetition pattern that let the
# hist_2/hist_3 relabeling typo slip into the LM pipeline unnoticed (see
# prepare_data_lm.r). theta0 values are kept as an explicit table (copied
# verbatim from the original script) rather than reconstructed
# algebraically, since unlike prepare_data_lm.r's order-swap logic there's
# no shared derivation to fall back on here -- these are just literal
# scenario definitions. Output file paths are unchanged from the original
# script.
normal_default_par <- list(n0 = c(100, 100, 100), n = 100, sigma2 = 1,
                            mu0 = 0, tau02 = 1, al = 1/2, bl = 1/2,
                            theta = 5, alpha = rep(1/4, 4))

normal_theta0 <- list(
  "1_1" = c(4.6, 4.8, 5),   "1_2" = c(3.6, 3.8, 4),   "1_3" = c(2.6, 2.8, 3),
  "2_1" = c(4.6, 5, 4.8),   "2_2" = c(3.6, 4, 3.8),   "2_3" = c(2.6, 3, 2.8),
  "3_1" = c(5, 4.8, 4.6),   "3_2" = c(4, 3.8, 3.6),   "3_3" = c(3, 2.8, 2.6),
  "4_1" = c(5, 5, 5),       "4_2" = c(4, 4, 4),       "4_3" = c(3, 3, 3),
  "5_1" = c(3, 7, 5),       "5_2" = c(5, 7, 3)
)

make_scenario <- function(theta0, par = normal_default_par) {
  modifyList(par, list(theta0 = theta0))
}

# run simulations
n_cores <- 15
n_sim <- 200

# Reproducibility: every scenario gets its own fixed seed (derived from its
# position in `normal_theta0`) so re-running this script reproduces the
# same 200 replicates per scenario. See the `seed` argument added to
# sim_sce() in aux_fun_sim.R -- pass seed = NULL here to go back to the old
# non-reproducible behaviour.
base_seed <- 20240521

for (i in seq_along(normal_theta0)) {
  key <- names(normal_theta0)[i]
  sce <- make_scenario(normal_theta0[[key]])
  file <- sprintf("results/samples/normal/sim%s_normal_fixed_var.qs2", key)
  message(sprintf("[%d/%d] sim_sce(normal_fixed_var): %s -> %s", i, length(normal_theta0), key, file))
  sim <- sim_sce(model = "normal_fixed_var", n_cores, n_sim, sce,
                 gamma_model_normal, delta_model_normal, seed = base_seed + i)
  qs2::qs_save(sim, file = file)
}
