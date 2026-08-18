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

# Compile the model
gamma_model_bin <- cmdstan_model("code/models/bin/gamma_bin.stan")
delta_model_bin <- cmdstan_model("code/models/bin/eta_bin.stan")

#-------------------------------------------------------------------------------
# Run simulations for the binomial model with different scenarios
#-------------------------------------------------------------------------------
#
# Every scenario shares the same n0/n/a/b/al/bl/alpha and differs only in
# theta0 (and, for the three "flat prior" groups, in which of a/b/al/bl/alpha
# are loosened). This used to be ~300 lines of copy-pasted
# `sceX.Y <- list(...)` / `simX.Y <- sim_sce(...)` / `qs_save(...)` blocks --
# the exact repetition pattern that let a scenario-relabeling typo slip into
# the LM pipeline unnoticed (see prepare_data_lm.r's "hist_1"/"hist3" bug).
# Scenarios are now defined once, in scenarios_bin.R, and run through a
# single loop below. That file is shared with analysis_bin.r, which needs
# the same scenarios' prior hyperparameters (scenarios[["3_1_all_flat"]]$a
# etc.) for a couple of closed-form prior plots -- keeping one source of
# truth means the two scripts can no longer drift out of sync. Output file
# paths are unchanged from the original script, so every downstream
# qs2::qs_read("results/samples/bin/...") call in analysis_bin.r /
# analysis_seq_bin.r / get_summary_bin.r still resolves.
source("code/sim_exp/bin/scenarios_bin.R")

# run simulations
n_cores <- 15
n_sim <- 200

# Reproducibility: every scenario gets its own fixed seed (derived from its
# position in `scenarios`) so re-running this script reproduces the same
# 200 replicates per scenario. See the `seed` argument added to sim_sce()
# in aux_fun_sim.R -- pass seed = NULL here to go back to the old
# non-reproducible behaviour.
base_seed <- 20240521
sims <- vector("list", length(scenarios))
names(sims) <- names(scenarios)

for (i in seq_along(scenarios)) {
  key <- names(scenarios)[i]
  message(sprintf("[%d/%d] sim_sce(bin): %s -> %s", i, length(scenarios), key, scenario_files[[key]]))
  sims[[key]] <- sim_sce(model = "bin", n_cores, n_sim, scenarios[[key]],
                          gamma_model_bin, delta_model_bin,
                          seed = base_seed)
  qs2::qs_save(sims[[key]], file = scenario_files[[key]])
}

#-------------------------------------------------------------------------------
# Prior-only draws (post = 0) for the three "all flat" scenario III specs
#-------------------------------------------------------------------------------

for (id in names(sce3_theta0)) {
  key <- paste0(id, "_all_flat")
  file <- paste0("results/samples/bin/samples_prior_sce", id, "_all_flat.qs2")
  message("sample_sce_bin (prior only): ", key, " -> ", file)
  set.seed(base_seed + which(names(scenarios) == key))
  samples_prior <- sample_sce_bin(c(scenarios[[key]], post = 0), gamma_model_bin, delta_model_bin)
  qs2::qs_save(samples_prior, file = file)
}
