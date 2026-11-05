# =============================================================================
# SAMPLING step for the NPP/ONPP/Logit-Normal prior comparison, linear-model
# (lm) family -- restricted to the "varying beta, n0 > n" scenario section
# only (12 scenarios: Order I-IV x Congruence I-III), matching
# analysis_seq_lm.r's sce_list/df_sce. The "varying sigma" and "n0 < n"
# sections that also exist in this family (see simulations_lm.r) are NOT
# covered here -- add them the same way if needed later.
#
# Fits NPP (i.i.d. Beta, eta_lm.stan), ONPP (order-restricted, gamma_lm.stan),
# and the hierarchical compound-symmetry Logit-Normal prior
# (logitnorm_lm.stan), across the 12 scenarios x n_rep replicates, and
# caches every draw to disk.
#
# Unlike the bin family (prior_comparison_sample.R), this script does NOT
# simulate data itself: the historical + current data for every replicate
# already exists on disk, prepared by prepare_data_lm.r
# (results/sim_data/lm/sce<Order>_<Congruence>_varying_beta_n0_ge_n.qs2,
# one file per scenario, each holding up to 200 replicates' worth of
# ready-to-use Stan data). This script just reads that data and fits the 3
# priors on top of it.
#
# This script does NO plotting and builds NO tables -- see
# prior_comparison_lm_plot.R for that (it reads the cache this script
# writes and does zero further Stan sampling). logitnorm_aux_lm_plot.R and
# logitnorm_aux_lm_scenario_plot.R also read this same cache, for the
# hierarchical prior's auxiliary variables (mu/tau/rho/z/e).
#
# Run from the project root. Delete
# results/samples/lm/prior_comparison_lm_replicates.qs2 to force a full
# re-fit (e.g. after changing the hyperparameters or n_rep below).
# =============================================================================

library(cmdstanr)
library(posterior)
library(parallel)
library(pbapply)
library(qs2)

#-------------------------------------------------------------------------------
# Compile models
#-------------------------------------------------------------------------------

npp_model       <- cmdstan_model("code/models/lm/eta_lm.stan")       # NPP / NPP-SEQ
onpp_model      <- cmdstan_model("code/models/lm/gamma_lm.stan")     # ONPP / ONPP-SEQ (reference)
logitnorm_model <- cmdstan_model("code/models/lm/logitnorm_lm.stan") # coupled, no order (hierarchical)

#-------------------------------------------------------------------------------
# Scenario grid: Order (Compatible/Semi-compat./Incompatible/Neutral) x
# Discrepancy (None/Small/Large), "varying beta, n0 > n" section only.
# Eps/Order/Discrepancy labels copied verbatim from analysis_seq_lm.r's
# df_sce so results here line up directly with its figures. Each scenario's
# data is read from disk (prepare_data_lm.r's output), not simulated.
#-------------------------------------------------------------------------------

lm_par <- list(p = 3, a = 2, b = 1, V0 = 100 * diag(3), mu0 = rep(0, 3))

df_sce <- data.frame(
  Eps = c(
    "(1.0, 0.5, 0.0)", "(1.5, 1.0, 0.5)", "(2.0, 1.5, 1.0)",
    "(1.0, 0.0, 0.5)", "(1.5, 0.5, 1.0)", "(2.5, 1.5, 2.0)",
    "(0.0, 0.5, 1.0)", "(0.5, 1.0, 1.5)", "(1.0, 1.5, 2.0)",
    "(0.0, 0.0, 0.0)", "(1.0, 1.0, 1.0)", "(2.0, 2.0, 2.0)"
  ),
  Order       = c(rep("Compatible",   3), rep("Semi-compat.", 3),
                   rep("Incompatible", 3), rep("Neutral",      3)),
  Discrepancy = rep(c("None", "Small", "Large"), 4)
)

sce_keys <- character(0)
sce_list <- list()
idx <- 1
for (order_i in 1:4) {
  for (cong_j in 1:3) {
    order_roman <- as.character(as.roman(order_i))
    cong_roman  <- as.character(as.roman(cong_j))
    key <- paste0(order_roman, "_", cong_roman)
    sce_file <- sprintf("results/sim_data/lm/sce%s_%s_varying_beta_n0_ge_n.qs2", order_roman, cong_roman)
    sce_list[[key]] <- list(
      file = sce_file,
      eps = df_sce$Eps[idx],
      order = df_sce$Order[idx],
      discrepancy = df_sce$Discrepancy[idx]
    )
    sce_keys <- c(sce_keys, key)
    idx <- idx + 1
  }
}

# Hyperparameters for the logistic-normal prior -- see logitnorm_lm.stan's
# header comment for what each knob controls. Same defaults as the bin
# family's logitnorm_hyper (prior_comparison_sample.R) -- the (0,1)-scale
# logit-normal behaviour these control doesn't depend on which likelihood
# family eta is feeding into.
logitnorm_hyper <- list(
  mu_mean = 0, mu_sd = 1.5,        # mu ~ Normal(0, 1.5) on the logit scale
  tau_sd  = 1,                     # tau ~ Normal+(0, 1) (half-normal)
  rho_a = 2, rho_b = 2             # rho ~ Beta(2,2): symmetric, avoids the 0/1 boundary
)

#-------------------------------------------------------------------------------
# Convert one replicate's prepared data (prepare_data_stan()'s per-replicate
# entry: X0_flat/dims_X0/startid_X0, y0_flat/dims_y0/startid_y0, X, y) into
# the Stan data fields shared by all 3 models.
#-------------------------------------------------------------------------------

base_stan_data <- function(rep_data, par) {
  list(
    K = nrow(rep_data$dims_X0),
    n0 = rep_data$dims_X0[, 1],
    n = nrow(rep_data$X),
    p = par$p,
    len_X0 = length(rep_data$X0_flat),
    len_y0 = length(rep_data$y0_flat),
    X0 = rep_data$X0_flat,
    startid_X0 = rep_data$startid_X0,
    y0 = rep_data$y0_flat,
    startid_y0 = rep_data$startid_y0,
    X = rep_data$X,
    y = as.vector(rep_data$y),
    a = par$a,
    b = par$b,
    V0 = par$V0,
    mu0 = par$mu0
  )
}

#-------------------------------------------------------------------------------
# Fit helper -- same chains/parallel_chains split as the bin family's
# prior_comparison_sample.R: chains run sequentially within a replicate,
# parallelism happens across (scenario, replicate) jobs instead.
#-------------------------------------------------------------------------------

fit_variant <- function(model, data, sampler_args = list()) {
  defaults <- list(chains = 4, parallel_chains = 1,
                    iter_warmup = 1000, iter_sampling = 1000,
                    adapt_delta = 0.99, max_treedepth = 12, refresh = 0)
  args <- utils::modifyList(defaults, sampler_args)
  do.call(model$sample, c(list(data = data), args))
}

#-------------------------------------------------------------------------------
# Fit all 3 priors (NPP, ONPP, Logit-Normal hierarchical) on one
# replicate's data, at a given `post` stage. post = 1: usual posterior
# draws. post = 0: prior-only draws (the replicate's data has no effect in
# that case -- see logitnorm_lm.stan's `if (post == 1)` branch -- kept in
# the data list only so the "prior" fit uses the same data-list shape as
# the "posterior" fit).
#
# `theta` (from generated quantities: beta_draw appended with sigma2_draw,
# a length p+1 vector) is kept whole rather than split into beta/sigma2
# here, matching analysis_seq_lm.r's draws$theta_npp[, 1:3] / [, 4]
# convention.
#-------------------------------------------------------------------------------

fit_replicate_draws <- function(rep_data, par, post = 1) {
  base <- base_stan_data(rep_data, par)

  data_npp       <- c(base, list(al = 1, bl = 1, post = post, seq = 0))
  data_onpp      <- c(base, list(alpha = rep(1, 4), post = post, seq = 0))
  data_logitnorm <- c(base, logitnorm_hyper, list(post = post, seq = 0))

  fit_npp   <- fit_variant(npp_model,       data_npp)
  fit_onpp  <- fit_variant(onpp_model,      data_onpp)
  fit_logit <- fit_variant(logitnorm_model, data_logitnorm)

  list(
    eta_npp   = as_draws_matrix(fit_npp$draws("eta")),
    eta_onpp  = as_draws_matrix(fit_onpp$draws("eta")),
    eta_logit = as_draws_matrix(fit_logit$draws("eta")),
    theta_npp   = as_draws_matrix(fit_npp$draws("theta")),
    theta_onpp  = as_draws_matrix(fit_onpp$draws("theta")),
    theta_logit = as_draws_matrix(fit_logit$draws("theta")),
    # Auxiliary variables of the hierarchical logit-normal prior (mu, tau,
    # rho, z, e[1:K]) -- see prior_comparison_sample.R's (bin family)
    # identical comment. Kept so logitnorm_aux_lm_plot.R / logitnorm_aux_
    # lm_scenario_plot.R can reuse this exact fit instead of re-fitting
    # logitnorm_lm.stan.
    aux_logit = as_draws_matrix(fit_logit$draws(c("mu", "tau", "rho", "z", "e"))),
    divergences = c(npp   = sum(fit_npp$sampler_diagnostics()[, , "divergent__"]),
                     onpp  = sum(fit_onpp$sampler_diagnostics()[, , "divergent__"]),
                     logit = sum(fit_logit$sampler_diagnostics()[, , "divergent__"]))
  )
}

fit_one_replicate <- function(rep_data, par) {
  fit_replicate_draws(rep_data, par, post = 1)
}

#-------------------------------------------------------------------------------
# Cache -- a single qs2 file with three parts: sce_list (scenario metadata,
# for downstream labelling), replicates_by_scenario (the main grid of
# post = 1 fits), and prior_cell (one extra post = 0 fit for a single
# chosen (scenario, replicate), used by the "prior vs. posterior" plot).
# Same structure as the bin family's cache (prior_comparison_sample.R).
#-------------------------------------------------------------------------------

cache_file <- "results/samples/lm/prior_comparison_lm_replicates.qs2"
cached <- if (file.exists(cache_file)) qs2::qs_read(cache_file) else list()

n_rep <- 100
n_cores <- max(1, parallel::detectCores() - 1)
base_seed <- 20261001

if (is.null(cached$replicates_by_scenario)) {
  # Read every scenario's prepared data up front (cheap -- these are just
  # qs2 reads, not Stan fits) so the cluster workers only need the already-
  # sliced-to-n_rep replicate lists, not the full on-disk files.
  sce_data_list <- lapply(sce_list, function(meta) {
    sce_raw <- qs2::qs_read(meta$file)
    n_avail <- length(sce_raw$data)
    list(par = sce_raw$par, reps = sce_raw$data[seq_len(min(n_rep, n_avail))])
  })

  cl <- makeCluster(n_cores)
  on.exit(stopCluster(cl), add = TRUE)

  clusterSetRNGStream(cl, base_seed)

  clusterEvalQ(cl, {
    library(cmdstanr)
    library(posterior)
  })

  clusterExport(cl,
                 varlist = c("npp_model", "onpp_model", "logitnorm_model",
                             "logitnorm_hyper",
                             "base_stan_data", "fit_variant", "fit_replicate_draws", "fit_one_replicate"),
                 envir = environment())

  job_grid <- do.call(rbind, lapply(names(sce_data_list), function(key) {
    data.frame(scenario = key, rep = seq_along(sce_data_list[[key]]$reps),
               stringsAsFactors = FALSE)
  }))

  pbapply::pboptions(type = "timer")
  raw_results <- pbapply::pblapply(
    seq_len(nrow(job_grid)),
    function(i, job_grid, sce_data_list) {
      tryCatch({
        key <- job_grid$scenario[i]
        r   <- job_grid$rep[i]
        fit_one_replicate(sce_data_list[[key]]$reps[[r]], sce_data_list[[key]]$par)
      }, error = function(e) list(.failed = TRUE, .error = conditionMessage(e)))
    },
    job_grid = job_grid, sce_data_list = sce_data_list,
    cl = cl
  )

  # Drop and report any failed jobs before combining results (mirrors
  # prior_comparison_sample.R's (bin family) failed-replicate handling).
  failed <- vapply(raw_results, function(x) isTRUE(x$.failed), logical(1))
  if (any(failed)) {
    msgs <- unique(vapply(raw_results[failed], `[[`, character(1), ".error"))
    warning(sprintf("%d/%d (scenario, replicate) jobs failed and were dropped:\n%s",
                     sum(failed), nrow(job_grid), paste(msgs, collapse = "\n")))
  }

  replicates_by_scenario <- lapply(names(sce_data_list), function(key) {
    raw_results[!failed & job_grid$scenario == key]
  })
  names(replicates_by_scenario) <- names(sce_data_list)

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
# One extra (scenario, replicate) cell fit at post = 0 (prior-only). Default
# cell: "III_I" (Incompatible, Discrepancy None) -- the lm-family analogue
# of the bin family's "Incompatible_None" default (the scenario class that
# breaks ONPP's order restriction -- see prior_comparison_lm_plot.R's
# header comment). Change prior_cell_scenario_key/prior_cell_rep_idx and
# re-run this script to fit a different cell.
#-------------------------------------------------------------------------------

prior_cell_scenario_key <- "III_I"
prior_cell_rep_idx <- 1
stopifnot(prior_cell_scenario_key %in% names(sce_list))
stopifnot(prior_cell_rep_idx >= 1,
          prior_cell_rep_idx <= length(replicates_by_scenario[[prior_cell_scenario_key]]))

need_prior_cell <- is.null(cached$prior_cell) ||
  !identical(cached$prior_cell$scenario_key, prior_cell_scenario_key) ||
  !identical(cached$prior_cell$rep_idx, prior_cell_rep_idx)

if (need_prior_cell) {
  sce_raw <- qs2::qs_read(sce_list[[prior_cell_scenario_key]]$file)
  rep_data <- sce_raw$data[[prior_cell_rep_idx]]

  message(sprintf("Fitting prior (post = 0) draws for %s, replicate #%d...",
                   prior_cell_scenario_key, prior_cell_rep_idx))
  prior_draws <- fit_replicate_draws(rep_data, sce_raw$par, post = 0)

  cached$prior_cell <- list(scenario_key = prior_cell_scenario_key,
                             rep_idx = prior_cell_rep_idx,
                             draws = prior_draws)
  qs2::qs_save(cached, file = cache_file)
} else {
  message("prior_cell already cached for ", prior_cell_scenario_key,
          " replicate #", prior_cell_rep_idx, " -- skipping refit.")
}

message("Done. Cache written to ", cache_file,
        " -- run prior_comparison_lm_plot.R (or logitnorm_aux_lm_plot.R /",
        " logitnorm_aux_lm_scenario_plot.R) to build figures/tables from it.")
