# =============================================================================
# SAMPLING step for the NPP/ONPP/HLNNPP (hierarchical logistic-normal)
# prior comparison on REAL data, mirroring code/sim_exp/bin/
# prior_comparison_sample.R but for two published clinical-trial datasets
# instead of simulated scenarios. Uses the project's own Stan models
# (eta_bin.stan / gamma_bin.stan / logitnorm_bin.stan) -- NOT the
# bespoke Metropolis-Hastings samplers in the user-supplied Script-Example.R
# -- so results aren't numerically identical to that script's printed
# output, but should be close for NPP/ONPP (same likelihood/prior, just a
# different sampler), since both target the same posterior.
#
# Datasets (both self-contained -- no external data file needed):
#   - "lupus": pediatric SLE trial, K = 2 historical datasets (two adult
#     lupus trials). n0 = c(275, 287), s0 = c(92, 125); current n = 39,
#     s = 17. Matches Script-Example.R's BerMNPP_MCMC*/BerOMNPP_MCMC*
#     calls for this example.
#   - "hovon": HOVON AML trial, K = 4 historical datasets. n0 = c(359, 252,
#     693, 437), s0 = c(279, 208, 598, 358); current n = 259, s = 214.
#
# The THIRD example in Script-Example.R (Docetaxel NSCLC trial with
# propensity-score stratification) needs patient-level data loaded from
# "D:/Bayesian/NPP-Order/code/NSCLC_simData.RData" -- a file that isn't
# part of this repo (confirmed: no NSCLC_simData.RData anywhere under the
# project folder) -- so it's skipped here. Drop that .RData file under
# this project (e.g. data/NSCLC_simData.RData) and extend `datasets` below
# with its n0/s0/n/s (summed per propensity-score stratum, as that script
# computes them) to add it.
#
# Hyperparameters:
#   - NPP (eta_bin.stan): al = bl = 0.5 -- i.i.d. Beta(1/2, 1/2) prior on
#     each eta[k], matching Script-Example.R's prior_delta_alpha/beta =
#     rep(1/2, K) for BerMNPP_MCMC1/2.
#   - ONPP (gamma_bin.stan): alpha = rep(1/(K+1), K+1) -- matches
#     Script-Example.R's prior_gamma (rep(1,3)/3 for lupus, rep(1/5,5) for
#     HOVON are exactly 1/(K+1) each).
#   - Current-data Beta(a, b) baseline: a = b = 0.5, matching
#     Script-Example.R's prior_p = c(1/2, 1/2) everywhere.
#   - HLNNPP (logitnorm_bin.stan): same default hyperprior used in
#     code/sim_exp/bin/prior_comparison_sample.R (mu_mean = 0, mu_sd = 1.5,
#     tau_sd = 1, rho_a = rho_b = 2) -- Script-Example.R has no HLNNPP
#     analogue to match, so this is the project's own standard choice.
#
# A "no-borrowing" reference posterior (current-data-only Beta(a+s, b+n-s))
# is also cached per dataset, using Script-Example.R's own reference
# numbers where it computed one (lupus_t/HOVON_t -- fit on the FULL trial
# sample sizes, not the interim n/s used for borrowing: 28/53 for lupus,
# 211/252 for HOVON -- see that script's comments near `lupus_t`/
# `HOVON_t`), so the "how much does borrowing move you toward the full
# trial's own answer" comparison is directly meaningful.
#
# Run from the project root. Delete results/samples/real_data/
# real_data_fits.qs2 to force a refit.
# =============================================================================

library(cmdstanr)
library(posterior)
library(qs2)

#-------------------------------------------------------------------------------
# Compile models (same 3 Stan files the simulated-scenario scripts use).
#-------------------------------------------------------------------------------

npp_model       <- cmdstan_model("code/models/bin/eta_bin.stan")       # NPP
onpp_model      <- cmdstan_model("code/models/bin/gamma_bin.stan")     # ONPP
logitnorm_model <- cmdstan_model("code/models/bin/logitnorm_bin.stan") # HLNNPP

logitnorm_hyper <- list(mu_mean = 0, mu_sd = 1.5, tau_sd = 1, rho_a = 2, rho_b = 2)

#-------------------------------------------------------------------------------
# Datasets -- n0/s0/n/s copied verbatim from Script-Example.R's calls to
# BerMNPP_MCMC1/BerOMNPP_MCMC1 for each example; ref_a/ref_b are that
# script's own "full trial, no borrowing" reference Beta parameters.
#-------------------------------------------------------------------------------

datasets <- list(
  lupus = list(
    label = "Pediatric SLE trial", K = 2,
    n0 = c(275, 287), s0 = c(92, 125), n = 39, s = 17,
    a = 0.5, b = 0.5, al = 0.5, bl = 0.5,
    ref_a = 28 + 0.5, ref_b = 53 - 28 + 0.5
  ),
  hovon = list(
    label = "HOVON AML trial", K = 4,
    n0 = c(359, 252, 693, 437), s0 = c(279, 208, 598, 358), n = 259, s = 214,
    a = 0.5, b = 0.5, al = 0.5, bl = 0.5,
    ref_a = 211 + 0.5, ref_b = 252 - 211 + 0.5
  )
)

#-------------------------------------------------------------------------------
# Fit helper -- same defaults as prior_comparison_sample.R's fit_variant().
#-------------------------------------------------------------------------------

fit_variant <- function(model, data, sampler_args = list()) {
  defaults <- list(chains = 4, parallel_chains = 4,
                    iter_warmup = 1000, iter_sampling = 1000,
                    adapt_delta = 0.99, max_treedepth = 12, refresh = 0)
  args <- utils::modifyList(defaults, sampler_args)
  do.call(model$sample, c(list(data = data), args))
}

#-------------------------------------------------------------------------------
# One dataset: fit all 3 priors' posterior (post = 1, seq = 0 -- the plain,
# non-sequential NPP/ONPP/HLNNPP variants the user asked to compare), plus
# the no-borrowing reference draws.
#-------------------------------------------------------------------------------

fit_dataset <- function(ds) {
  data_npp <- list(K = ds$K, n0 = ds$n0, s0 = ds$s0, n = ds$n, s = ds$s,
                    a = ds$a, b = ds$b, al = ds$al, bl = ds$bl,
                    post = 1, seq = 0)
  data_onpp <- list(K = ds$K, n0 = ds$n0, s0 = ds$s0, n = ds$n, s = ds$s,
                     a = ds$a, b = ds$b, alpha = rep(1 / (ds$K + 1), ds$K + 1),
                     post = 1, seq = 0)
  data_logitnorm <- c(list(K = ds$K, n0 = ds$n0, s0 = ds$s0, n = ds$n, s = ds$s,
                            a = ds$a, b = ds$b, post = 1, seq = 0),
                       logitnorm_hyper)

  message("  Fitting NPP...")
  fit_npp <- fit_variant(npp_model, data_npp)
  message("  Fitting ONPP...")
  fit_onpp <- fit_variant(onpp_model, data_onpp)
  message("  Fitting HLNNPP...")
  fit_logit <- fit_variant(logitnorm_model, data_logitnorm)

  nsample <- length(as.numeric(as_draws_matrix(fit_npp$draws("theta"))))

  list(
    eta_npp   = as_draws_matrix(fit_npp$draws("eta")),
    eta_onpp  = as_draws_matrix(fit_onpp$draws("eta")),
    eta_logit = as_draws_matrix(fit_logit$draws("eta")),
    theta_npp   = as.numeric(as_draws_matrix(fit_npp$draws("theta"))),
    theta_onpp  = as.numeric(as_draws_matrix(fit_onpp$draws("theta"))),
    theta_logit = as.numeric(as_draws_matrix(fit_logit$draws("theta"))),
    aux_logit = as_draws_matrix(fit_logit$draws(c("mu", "tau", "rho"))),
    theta_ref = rbeta(nsample, ds$ref_a, ds$ref_b),
    divergences = c(npp   = sum(fit_npp$sampler_diagnostics()[, , "divergent__"]),
                     onpp  = sum(fit_onpp$sampler_diagnostics()[, , "divergent__"]),
                     logit = sum(fit_logit$sampler_diagnostics()[, , "divergent__"]))
  )
}

#-------------------------------------------------------------------------------
# Fit (or load) and cache.
#-------------------------------------------------------------------------------

cache_file <- "results/samples/real_data/real_data_fits.qs2"

if (!file.exists(cache_file)) {
  fits <- list()
  for (key in names(datasets)) {
    message("Fitting dataset: ", key, " (", datasets[[key]]$label, ")...")
    fits[[key]] <- fit_dataset(datasets[[key]])
  }

  total_divergences <- Reduce(`+`, lapply(fits, `[[`, "divergences"))
  message("Total divergent transitions: ",
          paste(names(total_divergences), total_divergences, sep = "=", collapse = ", "))

  cached <- list(datasets = datasets, fits = fits)
  qs2::qs_save(cached, file = cache_file)
  message("Done. Cache written to ", cache_file,
          " -- run real_data_plot.R to build figures/tables from it.")
} else {
  message("Cache already exists at ", cache_file, " -- skipping refit. ",
          "Delete it to force a refit.")
}
