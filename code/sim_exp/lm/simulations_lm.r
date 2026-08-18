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
delta_model_lm <- cmdstan_model("code/models/lm/eta_lm.stan")

n_cores <- 15
n_sim <- 200

#-------------------------------------------------------------------------------
# Run simulations for the linear-model scenarios prepared by prepare_data_lm.r
#-------------------------------------------------------------------------------
#
# prepare_data_lm.r writes one sce<Order>_<Congruence>_<suffix>.qs2 file per
# (order I-IV) x (congruence I-III) combination, for each of four
# varying-parameter/n0-relation sections (suffix). This used to be ~300
# lines of copy-pasted `sceX.Y <- qs_read(...)` / `simX.Y <- sim_sce(...)` /
# `qs_save(...)` blocks, one per scenario, repeated once per section. Same
# fix as simulations_bin.R: one loop over the same order/congruence grid
# that prepare_data_lm.r uses, driven by the same roman-numeral (input) /
# arabic-numeral (output) naming convention, so the two files can't drift
# out of sync about how many scenarios exist. All input/output file paths
# are unchanged from the original script.
# Reproducibility: every scenario gets its own fixed seed (derived from its
# position in the order/congruence grid, offset by `seed_offset` so the four
# sections and the small-size block don't reuse the same seeds) so
# re-running this script reproduces the same 200 replicates per scenario.
# See the `seed` argument added to sim_sce() in aux_fun_sim.R -- pass
# seed = NULL to go back to the old non-reproducible behaviour.
base_seed <- 20240521

run_lm_scenarios <- function(suffix, n_cores, n_sim, seed_offset = 0) {
  for (order_i in 1:4) {
    for (cong_j in 1:3) {
      order_roman <- as.character(as.roman(order_i))
      cong_roman  <- as.character(as.roman(cong_j))
      sce_file <- sprintf("results/sim_data/lm/sce%s_%s_%s.qs2", order_roman, cong_roman, suffix)
      sim_file <- sprintf("results/samples/lm/sim%d_%d_lm_%s.qs2", order_i, cong_j, suffix)
      message(sprintf("simulations_lm: sce%s_%s (%s) -> %s", order_roman, cong_roman, suffix, sim_file))
      sce <- qs_read(sce_file)
      seed <- base_seed + seed_offset + (order_i - 1) * 3 + cong_j
      sim <- sim_sce(model = "lm", n_cores, n_sim, sce, gamma_model_lm, delta_model_lm, seed = seed)
      qs_save(sim, file = sim_file)
    }
  }
}

#------------------------------------------------------------------------------
# varying beta, n0 > n -- small size (Order I only: high vs no congruence)
#------------------------------------------------------------------------------

sce1.1 <- qs_read("results/sim_data/lm/sceI_I_varying_beta_n0_ge_n_small_size.qs2")
sce1.3 <- qs_read("results/sim_data/lm/sceI_III_varying_beta_n0_ge_n_small_size.qs2")
sim1.1 <- sim_sce(model = "lm", n_cores, n_sim, sce1.1, gamma_model_lm, delta_model_lm, seed = base_seed)
qs_save(sim1.1,
  file = "results/samples/lm/sim1_1_lm_varying_beta_n0_ge_n_small_size.qs2"
)
sim1.3 <- sim_sce(model = "lm", n_cores, n_sim, sce1.3, gamma_model_lm, delta_model_lm, seed = base_seed + 1)
qs_save(sim1.3,
  file = "results/samples/lm/sim1_3_lm_varying_beta_n0_ge_n_small_size.qs2"
)

#------------------------------------------------------------------------------
# varying beta, n0 > n
#------------------------------------------------------------------------------
run_lm_scenarios("varying_beta_n0_ge_n", n_cores, n_sim, seed_offset = 100)

#------------------------------------------------------------------------------
# varying beta, n0 < n
#------------------------------------------------------------------------------
run_lm_scenarios("varying_beta_n0_le_n", n_cores, n_sim, seed_offset = 200)

#------------------------------------------------------------------------------
# varying sigma, n0 > n
#------------------------------------------------------------------------------
run_lm_scenarios("varying_sigma_n0_ge_n", n_cores, n_sim, seed_offset = 300)

#------------------------------------------------------------------------------
# varying sigma, n0 < n
#------------------------------------------------------------------------------
run_lm_scenarios("varying_sigma_n0_le_n", n_cores, n_sim, seed_offset = 400)
