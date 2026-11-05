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
n_cores <- detectCores()-1
n_sim <- 1000

# Reproducibility: every scenario gets its own fixed seed (derived from its
# position in `scenarios`) so re-running this script reproduces the same
# 200 replicates per scenario. See the `seed` argument added to sim_sce()
# in aux_fun_sim.R -- pass seed = NULL here to go back to the old
# non-reproducible behaviour.
base_seed <- 20240521

#-------------------------------------------------------------------------------
# Shared-replicate groups
#-------------------------------------------------------------------------------
#
# Scenarios that share the same data-generating parameters -- n0, n, theta,
# and theta0 up to permutation -- can share their underlying Monte Carlo
# replicates (see build_shared_replicate_data() in aux_fun_sim.R). Stan-side
# prior overrides (a/b/al/bl/alpha, used by the "3_X_*_flat" variants) don't
# affect the data-generating process at all, so a flat-prior variant lands in
# the same group as its base counterpart whenever its theta0 matches -- e.g.
# "3_1" and "3_1_eta_flat"/"3_1_eta_gamma_flat"/"3_1_all_flat" all share
# theta0 = (0.50, 0.40, 0.30) and form one 4-member group.
#
# This generalizes the fix originally built just for the "cancellation
# effect" pair ("5_1"/"5_2", theta0 = (0.10,0.90,0.50) vs (0.50,0.90,0.10)):
# without shared data, NPP/NPP-SEQ's reported MSE/WIS differ between the two
# only by ordinary Monte Carlo noise, even though their eta prior
# (eta_bin.stan) is exchangeable across historical datasets and their *true*
# MSE/WIS must be identical -- see build_shared_replicate_data()'s header
# comment for the full rationale. ONPP/ONPP-SEQ (gamma_bin.stan, eta =
# cumulative_sum(gamma)) are genuinely order-sensitive by design and are NOT
# expected to match even with shared data.
#
# The resulting groups over the current scenario table:
#   {1_1, 2_1, 3_1, 3_1_eta_flat, 3_1_eta_gamma_flat, 3_1_all_flat}  theta0 = (.30,.40,.50)
#   {1_2, 2_2, 3_2, 3_2_eta_flat, 3_2_eta_gamma_flat, 3_2_all_flat}  theta0 = (.20,.30,.40)
#   {1_3, 2_3, 3_3, 3_3_eta_flat, 3_3_eta_gamma_flat, 3_3_all_flat}  theta0 = (.05,.15,.25)
#   {4_1}, {4_2}, {4_3}                                               singletons (theta0 all-equal)
#   {5_1, 5_2}                                                        theta0 = (.10,.50,.90)
group_key <- function(sce) {
  paste(paste(sce$n0, collapse = ","), sce$n, sce$theta,
        paste(sort(sce$theta0), collapse = ","), sep = "|")
}

keys <- names(scenarios)
groups <- split(keys, vapply(keys, function(k) group_key(scenarios[[k]]), character(1)))

sims <- vector("list", length(scenarios))
names(sims) <- names(scenarios)

for (grp_keys in groups) {
  # Which scenarios in this group still need to be run -- replaces the old
  # hand-maintained numeric `if (i < 19) next` skip, which had to be bumped
  # by hand after every batch and gave no protection against silently
  # re-skipping a scenario whose output file had been deleted or renamed.
  pending <- grp_keys[!file.exists(scenario_files[grp_keys])]
  if (length(pending) == 0) next

  if (length(pending) == 1) {
    # No sharing partner pending this pass -- run it the old way, with its
    # own independently-drawn data.
    key <- pending
    message(sprintf("[%s] sim_sce(bin): %s -> %s", key, key, scenario_files[[key]]))
    sims[[key]] <- sim_sce(model = "bin", n_cores, n_sim, scenarios[[key]],
                            gamma_model_bin, delta_model_bin,
                            seed = base_seed)
    qs2::qs_save(sims[[key]], file = scenario_files[[key]])
  } else {
    # >= 2 pending scenarios in this group share the same theta0 *set* --
    # build one shared pool of replicates (drawn by theta0 VALUE, not
    # position, and correctly relabeled per scenario) and reuse it across
    # all of them.
    message(sprintf("[shared] sim_sce(bin) group (%s): %s",
                     paste(pending, collapse = ", "),
                     paste(scenario_files[pending], collapse = ", ")))
    ref <- scenarios[[pending[1]]]
    shared <- build_shared_replicate_data(
      theta0_list = setNames(lapply(pending, function(k) scenarios[[k]]$theta0), pending),
      n0 = ref$n0, n = ref$n, theta = ref$theta,
      num_sim = n_sim, seed = base_seed
    )
    for (key in pending) {
      sims[[key]] <- sim_sce(model = "bin", n_cores, n_sim, scenarios[[key]],
                              gamma_model_bin, delta_model_bin,
                              replicate_data = shared[[key]])
      qs2::qs_save(sims[[key]], file = scenario_files[[key]])
    }
  }
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
