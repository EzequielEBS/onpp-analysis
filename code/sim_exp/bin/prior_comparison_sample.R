# =============================================================================
# SAMPLING step for the NPP/ONPP/Logit-Normal prior comparison across the
# full Order x Discrepancy scenario grid (Bernoulli experiments, Section
# 5.1). Fits NPP (i.i.d. Beta), ONPP (order-restricted), and the
# hierarchical compound-symmetry Logit-Normal prior (logitnorm_bin.stan),
# across 12 scenarios x n_rep replicates, and caches every draw to disk.
#
# This script does NO plotting and builds NO tables -- see
# prior_comparison_plot.R for that (it reads the cache this script writes
# and does zero further Stan sampling). logitnorm_aux_plot.R also reads this
# same cache, for the hierarchical prior's auxiliary variables (mu/tau/rho/
# z/e). Splitting it this way means the expensive part (thousands of Stan
# fits) only has to happen once, and any plotting/table tweak downstream is
# instant -- no re-fitting required.
#
# Run from the project root (same convention as simulations_bin.R). Delete
# results/samples/bin/prior_comparison_replicates.qs2 to force a full
# re-fit (e.g. after changing sce_list, the hyperparameters, or n_rep
# below).
# =============================================================================

library(cmdstanr)
library(posterior)
library(parallel)
library(pbapply)
library(qs2)

#-------------------------------------------------------------------------------
# Compile models
#-------------------------------------------------------------------------------

npp_model       <- cmdstan_model("code/models/bin/eta_bin.stan")       # NPP / NPP-SEQ
onpp_model      <- cmdstan_model("code/models/bin/gamma_bin.stan")     # ONPP / ONPP-SEQ (reference)
logitnorm_model <- cmdstan_model("code/models/bin/logitnorm_bin.stan") # coupled, no order (hierarchical)

#-------------------------------------------------------------------------------
# Scenario grid: Order (Compatible/Semi-compatible/Incompatible/Neutral) x
# Discrepancy (None/Small/Large). theta0 values are copied verbatim from
# scenarios_bin.R's main_theta0 (its "1_1".."4_3" entries) so results here
# line up directly with the manuscript's own Figure 3 / Section 5.1 numbers.
# order/discrepancy are stored on each scenario purely for labelling the
# plots downstream -- they don't affect simulation or fitting.
#-------------------------------------------------------------------------------

base_par <- list(n0 = c(100, 100, 100), n = 100,
                  a = 1, b = 1, al = 1, bl = 1,
                  alpha = rep(1, 4), theta = 0.50)

scenario_theta0 <- list(
  "Compatible_None"      = c(0.30, 0.40, 0.50),
  "Compatible_Small"     = c(0.20, 0.30, 0.40),
  "Compatible_Large"     = c(0.05, 0.15, 0.25),
  "Semi-compatible_None"  = c(0.30, 0.50, 0.40),
  "Semi-compatible_Small" = c(0.20, 0.40, 0.30),
  "Semi-compatible_Large" = c(0.05, 0.25, 0.15),
  "Incompatible_None"    = c(0.50, 0.40, 0.30),
  "Incompatible_Small"   = c(0.40, 0.30, 0.20),
  "Incompatible_Large"   = c(0.25, 0.15, 0.05),
  "Neutral_None"          = c(0.50, 0.50, 0.50),
  "Neutral_Small"         = c(0.30, 0.30, 0.30),
  "Neutral_Large"         = c(0.10, 0.10, 0.10)
)

sce_list <- lapply(names(scenario_theta0), function(key) {
  parts <- strsplit(key, "_")[[1]]
  modifyList(base_par, list(theta0 = scenario_theta0[[key]],
                             order = parts[1], discrepancy = parts[2]))
})
names(sce_list) <- names(scenario_theta0)

# Hyperparameters for the logistic-normal prior -- see the header comment in
# logitnorm_bin.stan for what each knob controls. Tune these the same way
# you already tune al/bl/alpha in scenarios_bin.R's flat_variants.
logitnorm_hyper <- list(
  mu_mean = 0, mu_sd = 1.5,        # mu ~ Normal(0, 1.5) on the logit scale
  tau_sd  = 1,                     # tau ~ Normal+(0, 1) (half-normal)
  rho_a = 2, rho_b = 2             # rho ~ Beta(2,2): symmetric, avoids the 0/1 boundary
)

K <- length(base_par$n0)

#-------------------------------------------------------------------------------
# Fit helper: one dataset, one model, "Default" (post = 1, seq = 0) variant
#
# parallel_chains = 1: chains run SEQUENTIALLY within a replicate. Parallelism
# happens one level up, across (scenario, replicate) jobs (see the PSOCK
# cluster below) -- same division of labour as sim_sce()/.fit_family() in
# aux_fun_sim.R, which always samples with cmdstanr's own (sequential)
# default and instead parallelizes the outer pblapply() over replicates.
#-------------------------------------------------------------------------------

fit_variant <- function(model, data, sampler_args = list()) {
  defaults <- list(chains = 4, parallel_chains = 1,
                    iter_warmup = 1000, iter_sampling = 1000,
                    adapt_delta = 0.99, max_treedepth = 12, refresh = 0)
  args <- utils::modifyList(defaults, sampler_args)
  do.call(model$sample, c(list(data = data), args))
}

#-------------------------------------------------------------------------------
# Fit all 3 priors (NPP, ONPP, Logit-Normal hierarchical) on one (z0, z)
# dataset, at a given `post` stage. post = 1: usual posterior draws.
# post = 0: prior-only draws (z0/z have no effect in that case -- see
# logitnorm_bin.stan's `if (post == 1)` branch -- kept in the data list
# only so the "prior" fit uses the same data-list shape as the "posterior"
# fit for the same replicate).
#-------------------------------------------------------------------------------

fit_replicate_draws <- function(sce, K, z0, z, post = 1) {
  data_npp <- list(K = K, n0 = sce$n0, s0 = z0, n = sce$n, s = z,
                    a = sce$a, b = sce$b, al = sce$al, bl = sce$bl,
                    post = post, seq = 0)
  data_onpp <- list(K = K, n0 = sce$n0, s0 = z0, n = sce$n, s = z,
                     a = sce$a, b = sce$b, alpha = sce$alpha,
                     post = post, seq = 0)
  data_logitnorm <- c(list(K = K, n0 = sce$n0, s0 = z0, n = sce$n, s = z,
                            a = sce$a, b = sce$b, post = post, seq = 0),
                       logitnorm_hyper)

  fit_npp   <- fit_variant(npp_model,       data_npp)
  fit_onpp  <- fit_variant(onpp_model,      data_onpp)
  fit_logit <- fit_variant(logitnorm_model, data_logitnorm)

  list(
    eta_npp   = as_draws_matrix(fit_npp$draws("eta")),
    eta_onpp  = as_draws_matrix(fit_onpp$draws("eta")),
    eta_logit = as_draws_matrix(fit_logit$draws("eta")),
    theta_npp   = as.numeric(as_draws_matrix(fit_npp$draws("theta"))),
    theta_onpp  = as.numeric(as_draws_matrix(fit_onpp$draws("theta"))),
    theta_logit = as.numeric(as_draws_matrix(fit_logit$draws("theta"))),
    # Auxiliary variables of the hierarchical logit-normal prior (mu, tau,
    # rho, z, e[1:K]). Kept so logitnorm_aux_plot.R's posterior
    # plots/table can reuse this exact fit instead of re-fitting
    # logitnorm_bin.stan.
    aux_logit = as_draws_matrix(fit_logit$draws(c("mu", "tau", "rho", "z", "e"))),
    divergences = c(npp   = sum(fit_npp$sampler_diagnostics()[, , "divergent__"]),
                     onpp  = sum(fit_onpp$sampler_diagnostics()[, , "divergent__"]),
                     logit = sum(fit_logit$sampler_diagnostics()[, , "divergent__"]))
  )
}

#-------------------------------------------------------------------------------
# One replicate: simulate data, fit all 3 priors' POSTERIOR (post = 1),
# return eta/theta/aux draws plus the (z0, z) used -- kept so a later
# "prior vs. posterior" fit for a chosen (scenario, replicate) cell can
# reuse the identical data list (see the prior_cell block below).
#-------------------------------------------------------------------------------

fit_one_replicate <- function(sce, K) {
  z0 <- vapply(1:K, function(k) rbinom(1, sce$n0[k], sce$theta0[k]), integer(1))
  z  <- rbinom(1, sce$n, sce$theta)
  draws <- fit_replicate_draws(sce, K, z0, z, post = 1)
  draws$z0 <- z0
  draws$z  <- z
  draws
}

#-------------------------------------------------------------------------------
# Cache -- a single qs2 file with three parts: sce_list (scenario metadata,
# for downstream labelling), replicates_by_scenario (the main grid of
# post = 1 fits), and prior_cell (one extra post = 0 fit for a single
# chosen (scenario, replicate), used by the "prior vs. posterior" plot).
# Each part is filled in independently and the cache re-saved after each,
# so re-running this script after the main grid is already cached only
# does the (cheap) prior_cell part, and changing prior_cell_scenario_key/
# prior_cell_rep_idx below and re-running only re-fits that one cell.
#-------------------------------------------------------------------------------

cache_file <- "results/samples/bin/prior_comparison_replicates.qs2"
cached <- if (file.exists(cache_file)) qs2::qs_read(cache_file) else list()

n_rep <- 100
n_cores <- max(1, parallel::detectCores() - 1)
base_seed <- 20260922

if (is.null(cached$replicates_by_scenario)) {
  # n_rep = 100: 12 scenarios x 100 replicates x 3 models = 3600 Stan fits
  # total, spread across n_cores workers, in parallel across a single PSOCK
  # cluster shared by all (scenario, replicate) jobs (one job's 4 Stan fits
  # per worker, chains sequential within each fit -- see fit_variant()'s
  # comment above). Same pattern as sim_sce() in aux_fun_sim.R:
  # clusterSetRNGStream() for reproducibility, the compiled models + helper
  # functions exported once, pbapply::pblapply() with a tryCatch() per job
  # so one bad fit doesn't lose the whole batch. The paper's own Bernoulli
  # experiments use n_sim = 1000; bump n_rep up to match if you want
  # publication-grade density (and delete the cache file to force a re-run
  # at the new n_rep).
  cl <- makeCluster(n_cores)
  on.exit(stopCluster(cl), add = TRUE)

  clusterSetRNGStream(cl, base_seed)

  clusterEvalQ(cl, {
    library(cmdstanr)
    library(posterior)
  })

  clusterExport(cl,
                 varlist = c("npp_model", "onpp_model", "logitnorm_model",
                             "sce_list", "K", "logitnorm_hyper",
                             "fit_variant", "fit_replicate_draws", "fit_one_replicate"),
                 envir = environment())

  job_grid <- expand.grid(scenario = names(sce_list), rep = seq_len(n_rep),
                           stringsAsFactors = FALSE)

  pbapply::pboptions(type = "timer")
  raw_results <- pbapply::pblapply(
    seq_len(nrow(job_grid)),
    function(i, job_grid, sce_list, K) {
      tryCatch(
        fit_one_replicate(sce_list[[job_grid$scenario[i]]], K),
        error = function(e) list(.failed = TRUE, .error = conditionMessage(e))
      )
    },
    job_grid = job_grid, sce_list = sce_list, K = K,
    cl = cl
  )

  # Drop and report any failed jobs before combining results (mirrors
  # sim_sce()'s failed-replicate handling in aux_fun_sim.R).
  failed <- vapply(raw_results, function(x) isTRUE(x$.failed), logical(1))
  if (any(failed)) {
    msgs <- unique(vapply(raw_results[failed], `[[`, character(1), ".error"))
    warning(sprintf("%d/%d (scenario, replicate) jobs failed and were dropped:\n%s",
                     sum(failed), nrow(job_grid), paste(msgs, collapse = "\n")))
  }

  # Reshape the flat job list back into replicates_by_scenario[[scenario]] =
  # list of that scenario's (surviving) replicate draws.
  replicates_by_scenario <- lapply(names(sce_list), function(key) {
    raw_results[!failed & job_grid$scenario == key]
  })
  names(replicates_by_scenario) <- names(sce_list)

  all_replicates <- unlist(replicates_by_scenario, recursive = FALSE)
  total_divergences <- Reduce(`+`, lapply(all_replicates, `[[`, "divergences"))
  message("Total divergent transitions across all (scenario, replicate) jobs: ",
          paste(names(total_divergences), total_divergences, sep = "=", collapse = ", "))
  if (any(total_divergences > 0)) {
    message("If a method shows many divergences, raise adapt_delta/max_treedepth ",
            "in fit_variant()'s defaults above, or tighten logitnorm_hyper ",
            "(e.g. smaller tau_sd, larger rho_a/rho_b).")
  }

  cached$sce_list <- sce_list
  cached$replicates_by_scenario <- replicates_by_scenario
  qs2::qs_save(cached, file = cache_file)
} else {
  sce_list <- cached$sce_list
  replicates_by_scenario <- cached$replicates_by_scenario
  message("Main scenario grid already cached (", length(replicates_by_scenario[[1]]),
          " replicates/scenario) -- skipping refit.")
}

#-------------------------------------------------------------------------------
# One extra (scenario, replicate) cell fit at post = 0 (prior-only), reusing
# that replicate's own (z0, z) data list (see fit_replicate_draws()'s
# header comment -- z0/z have no effect on the draws in that case, this
# just keeps the two fits' data lists identical in shape). Very fast (no
# likelihood term to evaluate). Change prior_cell_scenario_key/
# prior_cell_rep_idx and re-run this script to fit a different cell --
# e.g. "Incompatible_None" is the scenario that breaks ONPP's order
# restriction (see prior_comparison_sample.R's/prior_comparison_plot.R's
# header comments), which is why it's the default.
#-------------------------------------------------------------------------------

prior_cell_scenario_key <- "Incompatible_None"
prior_cell_rep_idx <- 1
stopifnot(prior_cell_scenario_key %in% names(sce_list))
stopifnot(prior_cell_rep_idx >= 1,
          prior_cell_rep_idx <= length(replicates_by_scenario[[prior_cell_scenario_key]]))

need_prior_cell <- is.null(cached$prior_cell) ||
  !identical(cached$prior_cell$scenario_key, prior_cell_scenario_key) ||
  !identical(cached$prior_cell$rep_idx, prior_cell_rep_idx)

if (need_prior_cell) {
  sce_chosen <- sce_list[[prior_cell_scenario_key]]
  chosen     <- replicates_by_scenario[[prior_cell_scenario_key]][[prior_cell_rep_idx]]

  message(sprintf("Fitting prior (post = 0) draws for %s, replicate #%d...",
                   prior_cell_scenario_key, prior_cell_rep_idx))
  prior_draws <- fit_replicate_draws(sce_chosen, K, chosen$z0, chosen$z, post = 0)

  cached$prior_cell <- list(scenario_key = prior_cell_scenario_key,
                             rep_idx = prior_cell_rep_idx,
                             draws = prior_draws)
  qs2::qs_save(cached, file = cache_file)
} else {
  message("prior_cell already cached for ", prior_cell_scenario_key,
          " replicate #", prior_cell_rep_idx, " -- skipping refit.")
}

message("Done. Cache written to ", cache_file,
        " -- run prior_comparison_plot.R (or logitnorm_aux_plot.R) to build",
        " figures/tables from it.")
