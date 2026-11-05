# =============================================================================
# Internal helpers shared by sample_sce_bin / sample_sce_poi /
# sample_sce_normal_fixed_var / sample_sce_lm.
#
# The four sample_sce_* functions below used to be ~90% duplicated: each
# built four near-identical Stan data lists (onpp/onppseq/npp/nppseq), ran
# four $sample() calls, pulled divergence counts through a slow
# as_draws_df() round trip, and assembled the same 18-field return list by
# hand. That duplication is also how the normal-model sd/variance bug and
# the Poisson-model "post" bug (both fixed below) went unnoticed for a
# while: four independent copies of the same logic can quietly drift apart.
#
# These helpers factor the shared plumbing out; each sample_sce_* function
# now only supplies what's actually family-specific: the Stan data fields
# and the posterior draw transform for theta.
#
# NOTE for maintainers: sim_sce() runs these on a PSOCK cluster and
# explicitly clusterExport()s the helper names below (see sim_sce()). If you
# add another shared helper here, add its name to that clusterExport() call
# too, or workers will fail with "could not find function ...".
# =============================================================================

# Build the four (onpp / onppseq / npp / nppseq) Stan data lists from a
# common base list plus the fields that differ between the "ordered" (alpha,
# Dirichlet-prior) and "non-ordered" (al/bl, Beta-prior) variants.
.build_variant_data <- function(base, alpha, al, bl) {
  list(
    onpp    = c(base, list(alpha = alpha, seq = 0)),
    onppseq = c(base, list(alpha = alpha, seq = 1)),
    npp     = c(base, list(al = al, bl = bl, seq = 0)),
    nppseq  = c(base, list(al = al, bl = bl, seq = 1))
  )
}

# Fit gamma_model (onpp/onppseq) and delta_model (npp/nppseq) on the four
# data variants and return the four cmdstanr fit objects.
#
# sampler_args can set chains / iter_warmup / iter_sampling / refresh, plus
# adapt_delta either as a single value (applied to all four fits) or as a
# named list c(onpp = , onppseq = , npp = , nppseq = ) for per-fit control
# (the normal-fixed-variance model needs a stricter adapt_delta for its
# onpp fit than the other three; see sample_sce_normal_fixed_var()). If
# adapt_delta is omitted entirely, cmdstanr's own default is used.
#
# Each $sample() call is wrapped in tryCatch() so a failure identifies which
# of the four fits broke, instead of an opaque error several call frames
# away. sim_sce() additionally makes a whole *replicate* skippable if one
# of its four fits fails, so a single bad replicate doesn't abort an entire
# batch of hundreds of simulations.
.fit_family <- function(gamma_model, delta_model,
                         data_onpp, data_onppseq, data_npp, data_nppseq,
                         sampler_args = list()) {
  defaults <- list(chains = 4, iter_warmup = 2000, iter_sampling = 2000, refresh = 0)
  adapt_delta <- sampler_args$adapt_delta
  base_args <- utils::modifyList(defaults, sampler_args[setdiff(names(sampler_args), "adapt_delta")])

  adapt_delta_for <- function(variant) {
    if (is.null(adapt_delta)) return(NULL)
    if (is.list(adapt_delta)) return(adapt_delta[[variant]])
    unname(adapt_delta)
  }

  run <- function(model, data, variant) {
    call_args <- c(list(data = data), base_args)
    ad <- adapt_delta_for(variant)
    if (!is.null(ad)) call_args$adapt_delta <- ad
    tryCatch(
      do.call(model$sample, call_args),
      error = function(e) {
        stop(sprintf("cmdstanr sampling failed for the '%s' fit: %s",
                      variant, conditionMessage(e)), call. = FALSE)
      }
    )
  }

  list(
    onpp    = run(gamma_model, data_onpp,    "onpp"),
    onppseq = run(gamma_model, data_onppseq, "onppseq"),
    npp     = run(delta_model, data_npp,     "npp"),
    nppseq  = run(delta_model, data_nppseq,  "nppseq")
  )
}

# Extract a named variable's draws (as a draws_matrix) for all four fits.
.extract_draws <- function(fits, variable) {
  lapply(fits, function(fit) fit$draws(variables = variable) %>% as_draws_matrix())
}

# Count divergent transitions directly from the sampler diagnostics array.
# (The original code round-tripped this through as_draws_df() just to sum
# one column -- summing the array slice directly is equivalent and faster.)
.count_divergences <- function(fit) {
  sum(fit$sampler_diagnostics()[, , "divergent__"])
}

# Assemble the standard sample_sce_* return value (hattheta/hatdelta means,
# divergence counts, and the raw draws) from the four fits plus the four
# eta/delta draws and four theta draws. theta_draws elements may be plain
# vectors (bin/poi/normal) or matrices (lm, one column per beta + sigma) --
# colMeans()/mean() is picked automatically per element.
.assemble_sce_result <- function(fits, eta_draws, theta_draws) {
  hattheta <- lapply(theta_draws, function(x) if (is.matrix(x)) colMeans(x) else mean(x))
  hatdelta <- lapply(eta_draws, colMeans)
  divergences <- lapply(fits, .count_divergences)

  list(
    hattheta_npp = hattheta$npp, hattheta_nppseq = hattheta$nppseq,
    hattheta_onpp = hattheta$onpp, hattheta_onppseq = hattheta$onppseq,
    hatdelta_npp = hatdelta$npp, hatdelta_nppseq = hatdelta$nppseq,
    hatdelta_onpp = hatdelta$onpp, hatdelta_onppseq = hatdelta$onppseq,
    divergences = list(divergences_npp = divergences$npp,
                       divergences_nppseq = divergences$nppseq,
                       divergences_onpp = divergences$onpp,
                       divergences_onppseq = divergences$onppseq),
    delta_npp = eta_draws$npp, delta_nppseq = eta_draws$nppseq,
    delta_onpp = eta_draws$onpp, delta_onppseq = eta_draws$onppseq,
    theta_npp = theta_draws$npp, theta_nppseq = theta_draws$nppseq,
    theta_onpp = theta_draws$onpp, theta_onppseq = theta_draws$onppseq
  )
}

# =============================================================================
# function to simulate binomial data and run the four models (NPP / NPP-SEQ
# / ONPP / ONPP-SEQ)
# =============================================================================
sample_sce_bin <- function(par_list, gamma_model, delta_model, sampler_args = list()) {
  # get parameters
  n0 <- par_list$n0
  n <- par_list$n
  a <- par_list$a
  b <- par_list$b
  al <- par_list$al
  bl <- par_list$bl
  theta0 <- par_list$theta0
  theta <- par_list$theta
  alpha <- par_list$alpha
  K <- length(n0)
  z0 <- par_list$z0
  z <- par_list$z
  if (is.null(z0)) {
    z0 <- unlist(lapply(1:K, function(i) rbinom(1, n0[i], theta0[i])))
  }
  if (is.null(z)) {
    z <- rbinom(1, n, theta)
  }
  post <- if (is.null(par_list$post)) 1 else par_list$post

  # define data lists
  base <- list(a = a, b = b, K = K, n = n, s = z, n0 = n0, s0 = z0, post = post)
  data <- .build_variant_data(base, alpha, al, bl)

  # sample from the models
  fits <- .fit_family(gamma_model, delta_model,
                       data$onpp, data$onppseq, data$npp, data$nppseq,
                       utils::modifyList(list(adapt_delta = 0.98), sampler_args))

  eta_draws <- .extract_draws(fits, "eta")

  # FIX: theta is drawn directly from each Stan fit's generated quantities
  # block, which already branches correctly on post (posterior vs.
  # prior-only) and seq -- see eta_bin.stan / gamma_bin.stan. This used to
  # be recomputed here in R with a formula that always folds in the current
  # data (z, n), i.e. the post = 1 posterior-update formula, so a post = 0
  # (prior-only) call silently got posterior-updated theta draws instead of
  # prior-only ones. Mirrors the pattern already used in
  # sample_sce_normal_fixed_var().
  theta_draws <- lapply(fits, function(fit) {
    fit$draws(variables = "theta") %>% as_draws_df() %>% pull(theta)
  })

  .assemble_sce_result(fits, eta_draws, theta_draws)
}

# =============================================================================
# function to simulate Poisson data and run the four models
# =============================================================================
sample_sce_poi <- function(par_list, gamma_model, delta_model, sampler_args = list()) {
  # get parameters
  t0 <- par_list$t0
  t <- par_list$t
  a <- par_list$a
  b <- par_list$b
  al <- par_list$al
  bl <- par_list$bl
  theta0 <- par_list$theta0
  theta <- par_list$theta
  alpha <- par_list$alpha
  K <- length(t0)
  z0 <- par_list$z0
  z <- par_list$z
  if (is.null(z0)) {
    z0 <- unlist(lapply(1:K, function(i) rpois(1, t0[i] * theta0[i])))
  }
  if (is.null(z)) {
    z <- rpois(1, t * theta)
  }
  # FIX: this used to hardcode post = 1 in every data list regardless of
  # par_list$post, unlike the bin/normal/lm versions -- a prior-only
  # (post = 0) call for the Poisson model was silently ignored and always
  # sampled from the posterior instead.
  post <- if (is.null(par_list$post)) 1 else par_list$post

  # define data lists
  base <- list(a = a, b = b, K = K, t = t, z = z, t0 = t0, z0 = z0, post = post)
  data <- .build_variant_data(base, alpha, al, bl)

  # sample from the models
  fits <- .fit_family(gamma_model, delta_model,
                       data$onpp, data$onppseq, data$npp, data$nppseq,
                       utils::modifyList(list(adapt_delta = 0.98), sampler_args))

  eta_draws <- .extract_draws(fits, "delta")

  theta_draws <- list(
    npp     = rgamma(nrow(eta_draws$npp),
                      z + a + eta_draws$npp %*% z0,
                      b + eta_draws$npp %*% t0 + t),
    nppseq  = rgamma(nrow(eta_draws$nppseq),
                      z + K * (a - 1) + 1 + eta_draws$nppseq %*% z0,
                      K * b + eta_draws$nppseq %*% t0 + t),
    onpp    = rgamma(nrow(eta_draws$onpp),
                      z + a + eta_draws$onpp %*% z0,
                      b + eta_draws$onpp %*% t0 + t),
    onppseq = rgamma(nrow(eta_draws$onppseq),
                      z + K * (a - 1) + 1 + eta_draws$onppseq %*% z0,
                      K * b + eta_draws$onppseq %*% t0 + t)
  )

  .assemble_sce_result(fits, eta_draws, theta_draws)
}

# =============================================================================
# function to simulate normal (fixed-variance) data and run the four models
# =============================================================================
sample_sce_normal_fixed_var <- function(par_list, gamma_model, delta_model, sampler_args = list()) {
  # get parameters
  n0 <- par_list$n0
  n <- par_list$n
  mu0 <- par_list$mu0
  tau02 <- par_list$tau02
  al <- par_list$al
  bl <- par_list$bl
  alpha <- par_list$alpha
  theta0 <- par_list$theta0
  theta <- par_list$theta
  sigma2 <- par_list$sigma2
  post <- if (is.null(par_list$post)) 1 else par_list$post

  K <- length(n0)
  y0 <- par_list$y0
  y <- par_list$y
  # FIX: sigma2 is a *variance* (it's also passed to the Stan model as
  # sigma2), but this used to be plugged straight into rnorm(sd = sigma2)
  # without sqrt(). That's invisible whenever sigma2 = 1 (every scenario
  # file currently used sigma2 = 1), but it would silently simulate data
  # with the wrong spread for any other variance.
  if (is.null(y0)) {
    y0 <- lapply(1:K, function(i)
      rnorm(n0[i], mean = theta0[i], sd = sqrt(sigma2)))
  }
  if (is.null(y)) {
    y <- rnorm(n, mean = theta, sd = sqrt(sigma2))
  }

  # define data lists
  base <- list(
    K = K,
    n0 = n0,
    n = n,
    mu0 = mu0,
    tau02 = tau02,
    sigma2 = sigma2,
    ybar0 = unlist(lapply(1:K, function(i) mean(y0[[i]]))),
    ybar = mean(y),
    post = post
  )
  data <- .build_variant_data(base, alpha, al, bl)

  # sample from the models
  #
  # NOTE: the onpp fit uses a stricter adapt_delta (7 nines) than the other
  # three (6 nines) here -- that asymmetry was already present in the
  # original code. It's kept as-is rather than silently unified, since
  # doing so would be a real (if small) change to sampler behaviour; if the
  # difference wasn't intentional, tighten/loosen it explicitly here.
  fits <- .fit_family(gamma_model, delta_model,
                       data$onpp, data$onppseq, data$npp, data$nppseq,
                       utils::modifyList(
                         list(adapt_delta = list(onpp = 0.9999999,
                                                  onppseq = 0.999999,
                                                  npp = 0.999999,
                                                  nppseq = 0.999999)),
                         sampler_args
                       ))

  eta_draws <- .extract_draws(fits, "eta")
  theta_draws <- lapply(fits, function(fit) {
    fit$draws(variables = "theta") %>% as_draws_df() %>% pull(theta)
  })

  .assemble_sce_result(fits, eta_draws, theta_draws)
}

# =============================================================================
# function to run the four models for the linear-regression case (data is
# prepared ahead of time by prepare_data_stan(), unlike the other three
# families which simulate their own data)
# =============================================================================
sample_sce_lm <- function(par_list, data, gamma_model, delta_model, sampler_args = list()) {
  X0 <- data$X0_flat
  startid_X0 <- data$startid_X0
  n0 <- as.vector(data$dims_y0[, 1])
  y0 <- data$y0_flat
  startid_y0 <- data$startid_y0
  K <- length(startid_X0)
  X <- data$X
  p <- ncol(X)
  n <- nrow(X)
  len_X0 <- length(X0)
  len_y0 <- length(y0)
  y <- as.vector(data$y)
  # get parameters
  mu0 <- par_list$mu0
  V0 <- par_list$V0
  a <- par_list$a
  b <- par_list$b
  al <- par_list$al
  bl <- par_list$bl
  alpha <- par_list$alpha
  post <- if (is.null(par_list$post)) 1 else par_list$post

  # define data lists
  base <- list(
    K = K, n0 = n0, n = n, p = p,
    len_X0 = len_X0, len_y0 = len_y0,
    X0 = X0, startid_X0 = startid_X0,
    y0 = y0, startid_y0 = startid_y0,
    X = X, y = y,
    a = a, b = b, V0 = V0, mu0 = mu0,
    post = post
  )
  data_variants <- .build_variant_data(base, alpha, al, bl)

  # sample from the models -- no adapt_delta override here (matches the
  # original, which left it commented out), so cmdstanr's own default is
  # used unless the caller supplies one via sampler_args.
  fits <- .fit_family(gamma_model, delta_model,
                       data_variants$onpp, data_variants$onppseq,
                       data_variants$npp, data_variants$nppseq,
                       sampler_args)

  eta_draws <- .extract_draws(fits, "eta")
  theta_draws <- lapply(fits, function(fit) fit$draws(variables = "theta") %>% as_draws_matrix())

  .assemble_sce_result(fits, eta_draws, theta_draws)
}

# function to run simulations in parallel
sim_sce <- function(model, num_cores, num_sim, sce, gamma_model, delta_model, seed = NULL,
                    replicate_data = NULL) {
  # replicate_data: optional list of length num_sim (bin/poi/normal_fixed_var
  # only -- lm already threads its own per-replicate `data` through sce$data)
  # overriding par_list fields per replicate, e.g. replicate_data[[i]] =
  # list(z0 = c(...), z = ...). Use this via build_shared_replicate_data()
  # to force two scenarios that share the same *set* of theta0 values (just
  # assigned to different historical-dataset positions) to draw from
  # exactly the same underlying data per replicate, instead of each
  # independently simulating its own -- see build_shared_replicate_data()'s
  # header comment for why that matters and what it does/doesn't buy you.
  # Start a cluster with the desired number of cores. on.exit() guarantees
  # the cluster is torn down even if something below errors (previously a
  # failed pblapply() run left worker processes running).
  cl <- makeCluster(num_cores)
  on.exit(stopCluster(cl), add = TRUE)

  # Optional reproducibility: nothing in this codebase previously called
  # set.seed() anywhere, so re-running a simulation batch never gave the
  # same numbers twice. A plain set.seed() in the calling script wouldn't
  # be enough here since the actual random draws (rbinom/rpois/rnorm inside
  # sample_sce_*(), plus cmdstan's own sampling) happen on the *worker*
  # processes, which each start with their own independent RNG state.
  # clusterSetRNGStream() switches every worker to the L'Ecuyer-CMRG
  # generator and seeds it deterministically from `seed`, so the same
  # `seed` + `num_cores` + `num_sim` reproduces the same batch of
  # replicates. Leave `seed = NULL` (the default) to keep the old
  # non-reproducible behaviour.
  if (!is.null(seed)) {
    clusterSetRNGStream(cl, seed)
  }

  # Load only the packages sample_sce_*() actually needs on the workers.
  # (The original loaded the full tidyverse on every worker for what's
  # really just dplyr's %>%/pull(); tidyverse also pulls in ggplot2, tidyr,
  # readr, purrr, stringr and forcats, none of which are used here, and
  # that's meaningfully slower to load on every one of num_cores workers,
  # every time sim_sce() is called.)
  clusterEvalQ(cl, {
    library(cmdstanr)
    library(dplyr)
    library(posterior)
  })
  # Export the models to the cluster
  clusterExport(cl,
                varlist = c("gamma_model",
                            "delta_model"),
                envir = environment()
  )
  # Export the internal helpers sample_sce_*() relies on. These live in the
  # global environment (this file is source()'d at top level), not in
  # sim_sce()'s own call frame, so they need their own clusterExport() with
  # envir = globalenv(). If you add another shared helper near the top of
  # this file, add its name here too.
  clusterExport(cl,
                varlist = c(".build_variant_data", ".fit_family",
                            ".extract_draws", ".assemble_sce_result",
                            ".count_divergences"),
                envir = globalenv()
  )

  # Run simulations in parallel
  if (model == "bin"){
    sample_sce <- sample_sce_bin
  } else if (model == "poi"){
    sample_sce <- sample_sce_poi
  } else if (model == "normal_fixed_var"){
    sample_sce <- sample_sce_normal_fixed_var
  } else if (model == "lm"){
    sample_sce <- sample_sce_lm
  } else {
    stop("Invalid model type. Choose from 'bin', 'poi', 'normal_fixed_var', or 'lm'.")
  }

  pbapply::pboptions(type = "timer")
  # Each replicate is wrapped in tryCatch() so that one bad replicate (e.g.
  # a cmdstan sampling failure) doesn't abort the whole batch and lose every
  # already-completed replicate -- which, at n_sim in the hundreds and
  # several minutes per replicate, used to mean losing hours of compute for
  # one flaky fit. Failed replicates are dropped (with a warning) before
  # results are combined below; previously they used replicate(num_sim, sce,
  # simplify = FALSE) to build a list of num_sim identical copies of `sce`
  # just so pblapply() had something the right length to iterate over --
  # iterating over seq_len(num_sim) directly avoids that up-front copying
  # (the randomness comes from inside sample_sce(), not from `sce` itself).
  if (model == "lm") {
    rep_sce <- sce$data
    results <- pbapply::pblapply(
      rep_sce,
      function(data, par_list, gamma_model, delta_model) {
        tryCatch(
          sample_sce(par_list, data, gamma_model, delta_model),
          error = function(e) list(.failed = TRUE, .error = conditionMessage(e))
        )
      },
      par_list = sce$par,
      gamma_model = gamma_model,
      delta_model = delta_model,
      cl = cl
    )
  } else {
    if (!is.null(replicate_data) && length(replicate_data) != num_sim) {
      stop(sprintf("replicate_data has length %d but num_sim = %d -- need one entry per replicate.",
                    length(replicate_data), num_sim))
    }
    results <- pbapply::pblapply(
      seq_len(num_sim),
      function(i, sce, gamma_model, delta_model, replicate_data) {
        # When replicate_data[[i]] is supplied (e.g. a pre-generated z0/z
        # pair from build_shared_replicate_data()), it overrides those
        # fields of sce for this replicate only; sample_sce_*() already
        # skips drawing its own z0/z whenever they're non-NULL. Without
        # replicate_data this is just `sce`, unchanged from before.
        par_list_i <- if (is.null(replicate_data)) sce else utils::modifyList(sce, replicate_data[[i]])
        tryCatch(
          sample_sce(par_list_i, gamma_model = gamma_model, delta_model = delta_model),
          error = function(e) list(.failed = TRUE, .error = conditionMessage(e))
        )
      },
      sce = sce,
      gamma_model = gamma_model,
      delta_model = delta_model,
      replicate_data = replicate_data,
      cl = cl
    )
  }

  # Drop and report any failed replicates before combining results.
  failed <- vapply(results, function(x) isTRUE(x$.failed), logical(1))
  if (any(failed)) {
    msgs <- unique(vapply(results[failed], `[[`, character(1), ".error"))
    warning(sprintf("sim_sce(): %d/%d replicate(s) failed and were dropped:\n  - %s",
                     sum(failed), length(results), paste(msgs, collapse = "\n  - ")),
            call. = FALSE)
    results <- results[!failed]
  }
  if (length(results) == 0) {
    stop("sim_sce(): every replicate failed; nothing to combine.", call. = FALSE)
  }

  # combine results
  if (model == "lm") {
    hattheta <- list()
    for (i in 1:length(results[[1]]$hattheta_npp)){
      hattheta[[i]] <- data.frame(
        hattheta_npp = sapply(results, function(x) x$hattheta_npp[i]),
        hattheta_nppseq = sapply(results, function(x) x$hattheta_nppseq[i]),
        hattheta_onpp = sapply(results, function(x) x$hattheta_onpp[i]),
        hattheta_onppseq = sapply(results, function(x) x$hattheta_onppseq[i])
      )
    }
  } else{
    hattheta <- data.frame(
      hattheta_npp = sapply(results, function(x) x$hattheta_npp),
      hattheta_nppseq = sapply(results, function(x) x$hattheta_nppseq),
      hattheta_onpp = sapply(results, function(x) x$hattheta_onpp),
      hattheta_onppseq = sapply(results, function(x) x$hattheta_onppseq)
    )
  }


  divergences <- data.frame(divergences_npp = sum(sapply(results, function(x) x$divergences$divergences_npp)),
                            divergences_nppseq = sum(sapply(results, function(x) x$divergences$divergences_nppseq)),
                            divergences_onpp = sum(sapply(results, function(x) x$divergences$divergences_onpp)),
                            divergences_onppseq = sum(sapply(results, function(x) x$divergences$divergences_onppseq))
  )

  hatdelta <- list()
  for (i in 1:length(results[[1]]$hatdelta_npp)){
    hatdelta[[i]] <- data.frame(
      hatdelta_npp = sapply(results, function(x) x$hatdelta_npp[i]),
      hatdelta_nppseq = sapply(results, function(x) x$hatdelta_nppseq[i]),
      hatdelta_onpp = sapply(results, function(x) x$hatdelta_onpp[i]),
      hatdelta_onppseq = sapply(results, function(x) x$hatdelta_onppseq[i])
    )
  }

  return(list(hattheta = hattheta,
              hatdelta = hatdelta,
              divergences = divergences,
              delta = lapply(results, function(x) list(delta_npp = x$delta_npp,
                                                        delta_nppseq = x$delta_nppseq,
                                                        delta_onpp = x$delta_onpp,
                                                        delta_onppseq = x$delta_onppseq)
              ),
              theta = lapply(results, function(x) list(theta_npp = x$theta_npp,
                                                        theta_nppseq = x$theta_nppseq,
                                                        theta_onpp = x$theta_onpp,
                                                        theta_onppseq = x$theta_onppseq)
              )
  ))
}

# get_hat_par() was removed here: it was never called anywhere in this
# repository (confirmed by search) and referenced fields
# (draws_beta_*/draws_sigma_*) that none of the sample_sce_*() functions
# above actually produce -- calling it would have errored. It looks like
# leftover code from an earlier version of the pipeline. Recoverable from
# git history if it's still needed.

# function to plot boxplots
plot_boxplot <- function(draws, variable) {
  # create data frame (all four columns/elements always have the same
  # number of draws, so one length computation suffices)
  n_draws <- length(draws[,1])

  c <- data.frame(x=c(rep("NPP1",n_draws),
                      rep("NPP2",n_draws),
                      rep("ONPP1",n_draws),
                      rep("ONPP2",n_draws)),
                  y=c(draws[,1],
                      draws[,2],
                      draws[,3],
                      draws[,4])
  )

  c$x <- factor(c$x,levels=c('NPP1','NPP2','ONPP1','ONPP2'),ordered=TRUE)

  # plot
  box <- ggplot(c, aes(x=x, y=y, fill=x)) +
    theme_bw()+labs(x = NULL) +
    labs(y = variable)+
    stat_boxplot(geom ='errorbar', coef = 0.5, width = 0.25)+
    geom_boxplot(outlier.colour="black",
                 width = 0.4, outlier.size = 0.3, coef = 0.5)+
    stat_summary(fun = mean, colour="yellow3", geom="point",
                 shape=18, size=3, show.legend = FALSE)+
    theme(legend.position="none",
          axis.title.y=element_text(angle= -270,  face='bold', size=11))+
    scale_fill_manual(values=c("#999999","red2", "green4", "dodgerblue1"))+
    scale_x_discrete(labels=c(expression(NPP), expression(NPP[seq]),
                              expression(ONPP), expression(ONPP[seq]))
    )
  # (previously ended on the bare assignment `box <- ggplot(...)`, which only
  # "worked" as a return value because of R's last-expression-return rule --
  # an explicit return() is less of a landmine for the next person who
  # appends a line after this.)
  return(box)
}

# function to compute BCI
compute_bci <- function(list_draws, alpha){
  bci_npp <- lapply(list_draws, function(x) quantile(x$theta_npp, 
                                                     probs = c((1-alpha)/2, (1+alpha)/2)))
  bci_nppseq <- lapply(list_draws, function(x) quantile(x$theta_nppseq, 
                                                        probs = c((1-alpha)/2, (1+alpha)/2)))
  bci_onpp <- lapply(list_draws, function(x) quantile(x$theta_onpp,
                                                      probs = c((1-alpha)/2, (1+alpha)/2)))
  bci_onppseq <- lapply(list_draws, function(x) quantile(x$theta_onppseq,
                                                         probs = c((1-alpha)/2, (1+alpha)/2)))
  return(list(bci_npp = bci_npp,
              bci_nppseq = bci_nppseq,
              bci_onpp = bci_onpp,
              bci_onppseq = bci_onppseq
  ))
}

plot_sce_bin <- function(j,
                          sim_sces = get("sim_sces", envir = .GlobalEnv),
                          true_value = get("true_value", envir = .GlobalEnv)) {
  # sim_sces/true_value default to looking up globals of the same name, so
  # existing call sites like plot_sce_bin(j) keep working unchanged -- but
  # the dependency is now visible in the signature and can be overridden
  # (e.g. for testing) instead of being an invisible requirement on the
  # caller's global environment.
  sim <- sim_sces[[j]]
  sce <- as.roman(ceiling(j/3))
  sce <- ifelse(j %% 3 == 1, paste0(sce,".I"), ifelse(j %% 3 == 2, paste0(sce,".II"), paste0(sce,".III")))
  plots_delta <- lapply(seq_along(sim$hatdelta), function(i) {
    plot <- plot_boxplot(sim$hatdelta[[i]], bquote(hat(delta)[.(i)])) +
      coord_flip() +
      ggtitle(paste0("Scenario ", sce))

    if (i>1) {
      plot <- plot + scale_x_discrete(labels = NULL) + ggtitle(NULL)
    }
    return(plot)
  })
  plot_theta <- plot_boxplot(sim$hattheta, expression(hat(theta))) +
    geom_hline(yintercept = true_value, linetype = "dotted", color = "red", size = 1) +
    coord_flip()
  # Combine plots
  plot <- Reduce(`+`, plots_delta) + plot_theta + plot_layout(ncol = length(plots_delta) + 1) +
    scale_x_discrete(labels = NULL)
  return(plot)
}

# function to compute coverage
compute_coverage <- function(bcis, true_theta) {
  coverage_npp <- sum(sapply(bcis$bci_npp, function(x) x[1] <= true_theta & true_theta <= x[2]))/
    length(bcis$bci_npp)
  coverage_nppseq <- sum(sapply(bcis$bci_nppseq, function(x) x[1] <= true_theta & true_theta <= x[2]))/
    length(bcis$bci_nppseq)
  coverage_onpp <- sum(sapply(bcis$bci_onpp, function(x) x[1] <= true_theta & true_theta <= x[2]))/
    length(bcis$bci_onpp)
  coverage_onppseq <- sum(sapply(bcis$bci_onppseq, function(x) x[1] <= true_theta & true_theta <= x[2]))/
    length(bcis$bci_onppseq)
  return(list(coverage_npp = coverage_npp,
              coverage_nppseq = coverage_nppseq,
              coverage_onpp = coverage_onpp,
              coverage_onppseq = coverage_onppseq
  ))
}

compute_avg_len <- function(bcis) {
  bcis_npp <- do.call(rbind,bcis$bci_npp)
  bcis_nppseq <- do.call(rbind,bcis$bci_nppseq)
  bcis_onpp <- do.call(rbind,bcis$bci_onpp)
  bcis_onppseq <- do.call(rbind,bcis$bci_onppseq)
  return(list(
    npp = mean(bcis_npp[,2] - bcis_npp[,1]),
    nppseq = mean(bcis_nppseq[,2] - bcis_nppseq[,1]),
    onpp = mean(bcis_onpp[,2] - bcis_onpp[,1]),
    onppseq = mean(bcis_onppseq[,2] - bcis_onppseq[,1])
  ))
}

compute_wis <- function(list_draws, quantile_level, true_value) {
  N <- length(list_draws)
  q_npp <- do.call(rbind,
                   lapply(list_draws, function(x) quantile(x$theta_npp, 
                                                           probs = quantile_level))
  )
  q_nppseq <- do.call(rbind,
                      lapply(list_draws, function(x) quantile(x$theta_nppseq, 
                                                              probs = quantile_level))
  )
  q_onpp <- do.call(rbind,
                    lapply(list_draws, function(x) quantile(x$theta_onpp,
                                                            probs = quantile_level))
  )
  q_onppseq <- do.call(rbind,
                       lapply(list_draws, function(x) quantile(x$theta_onppseq,
                                                               probs = quantile_level))
  )
  wis_npp <- wis(rep(true_value, N), q_npp, quantile_level)
  wis_nppseq <- wis(rep(true_value, N), q_nppseq, quantile_level)
  wis_onpp <- wis(rep(true_value, N), q_onpp, quantile_level)
  wis_onppseq <- wis(rep(true_value, N), q_onppseq, quantile_level)
  return(list(hat_wis_npp = mean(wis_npp),
              hat_wis_nppseq = mean(wis_nppseq),
              hat_wis_onpp = mean(wis_onpp),
              hat_wis_onppseq = mean(wis_onppseq)
  ))
}

compute_wis_lm <- function(list_draws, quantile_level, true_value) {
  N <- length(list_draws$theta)
  q_npp <- lapply(1:ncol(list_draws$theta[[1]]$theta_npp), function(j) {
    do.call(rbind,
            lapply(list_draws$theta, function(x) {
              quantile(x$theta_npp[, j], probs = quantile_level)
            })
      )
  })
  q_nppseq <- lapply(1:ncol(list_draws$theta[[1]]$theta_nppseq), function(j) {
    do.call(rbind,
            lapply(list_draws$theta, function(x) {
              quantile(x$theta_nppseq[, j], probs = quantile_level)
            })
    )
  })
  q_onpp <- lapply(1:ncol(list_draws$theta[[1]]$theta_onpp), function(j) {
    do.call(rbind,
            lapply(list_draws$theta, function(x) {
              quantile(x$theta_onpp[, j], probs = quantile_level)
            })
    )
  })
  q_onppseq <- lapply(1:ncol(list_draws$theta[[1]]$theta_onppseq), function(j) {
    do.call(rbind,
            lapply(list_draws$theta, function(x) {
              quantile(x$theta_onppseq[, j], probs = quantile_level)
            })
    )
  })
  wis_npp <- unlist(lapply(1:length(q_npp), function(i) {
    mean(wis(rep(true_value[i], N), q_npp[[i]], quantile_level))
  }))
  wis_nppseq <- unlist(lapply(1:length(q_nppseq), function(i) {
    mean(wis(rep(true_value[i], N), q_nppseq[[i]], quantile_level))
  }))
  wis_onpp <- unlist(lapply(1:length(q_onpp), function(i) {
    mean(wis(rep(true_value[i], N), q_onpp[[i]], quantile_level))
  }))
  wis_onppseq <- unlist(lapply(1:length(q_onppseq), function(i) {
    mean(wis(rep(true_value[i], N), q_onppseq[[i]], quantile_level))
  }))
  return(list(hat_wis_npp = wis_npp,
              hat_wis_nppseq = wis_nppseq,
              hat_wis_onpp = wis_onpp,
              hat_wis_onppseq = wis_onppseq
  ))
}

# function to plot BCI
plot_bci <- function(bci, sim, true_val, prob, coverage, scenario, model, y_range = c(0, 1)) {
  # FIX: this used to hardcode the reference line at 0.5 and the y-axis
  # window at [0, 1] regardless of what model/scenario was being plotted --
  # harmless for the binomial case (theta is a probability, true_val is
  # already 0.5), but silently wrong for any model whose parameter doesn't
  # live in [0, 1] (e.g. the normal-model theta ~ 5): every point would be
  # clipped outside the fixed window and the dashed reference line wouldn't
  # correspond to the true value at all. The reference line now uses the
  # `true_val` argument the function already receives instead of a
  # duplicated literal, and the y-window is a parameter (`y_range`) that
  # defaults to the original c(0, 1) -- so every existing bin/poi call site
  # is unaffected -- with y_range = NULL letting ggplot auto-scale for
  # models where no fixed window makes sense.
  scale_y_col1 <- if (is.null(y_range)) NULL else scale_y_continuous(limits = y_range)
  scale_y_col2 <- if (is.null(y_range)) {
    scale_y_continuous(breaks = NULL)
  } else {
    scale_y_continuous(limits = y_range, breaks = NULL)
  }

  bci_npp <- data.frame(sample = 1:length(bci$bci_npp),
                        hattheta = sim$hattheta$hattheta_npp,
                        lower = unlist(lapply(bci$bci_npp, 
                                              function(x) x[[1]])
                        ),
                        upper = unlist(lapply(bci$bci_npp, 
                                              function(x) x[[2]])
                        ),
                        include = unlist(lapply(bci$bci_npp, 
                                                function(x) x[[1]] <= true_val & x[[2]] >= true_val)
                        )
  )
  bci_nppseq <- data.frame(sample = 1:length(bci$bci_nppseq),
                           hattheta = sim$hattheta$hattheta_nppseq,
                           lower = unlist(lapply(bci$bci_nppseq, 
                                                 function(x) x[[1]])
                           ),
                           upper = unlist(lapply(bci$bci_nppseq, 
                                                 function(x) x[[2]])
                           ),
                           include = unlist(lapply(bci$bci_nppseq, 
                                                   function(x) x[[1]] <= true_val & x[[2]] >= true_val)
                           )
  )
  bci_onpp <- data.frame(sample = 1:length(bci$bci_onpp),
                         hattheta = sim$hattheta$hattheta_onpp,
                         lower = unlist(lapply(bci$bci_onpp, 
                                               function(x) x[[1]])
                         ),
                         upper = unlist(lapply(bci$bci_onpp, 
                                               function(x) x[[2]])
                         ),
                         include = unlist(lapply(bci$bci_onpp, 
                                                 function(x) x[[1]] <= true_val & x[[2]] >= true_val)
                         )
  )
  bci_onppseq <- data.frame(sample = 1:length(bci$bci_onppseq),
                            hattheta = sim$hattheta$hattheta_onppseq,
                            lower = unlist(lapply(bci$bci_onppseq, 
                                                  function(x) x[[1]])
                            ),
                            upper = unlist(lapply(bci$bci_onppseq, 
                                                  function(x) x[[2]])
                            ),
                            include = unlist(lapply(bci$bci_onppseq, 
                                                    function(x) x[[1]] <= true_val & x[[2]] >= true_val)
                            )
  )
  
  bci_npp$include <- factor(bci_npp$include, 
                               levels = c(TRUE, FALSE))
  bci_nppseq$include <- factor(bci_nppseq$include, 
                                  levels = c(TRUE, FALSE))
  bci_onpp$include <- factor(bci_onpp$include,
                                levels = c(TRUE, FALSE))
  bci_onppseq$include <- factor(bci_onppseq$include,
                                   levels = c(TRUE, FALSE))
  
  # add dummy row
  bci_npp <- rbind(bci_npp, 
                   data.frame(sample = NA,
                              hattheta = NA,
                              lower = NA,
                              upper = NA,
                              include = factor("FALSE", levels = c("TRUE", "FALSE"))
                   )
  )
  bci_nppseq <- rbind(bci_nppseq, 
                    data.frame(sample = NA,
                               hattheta = NA,
                               lower = NA,
                               upper = NA,
                               include = factor("FALSE", levels = c("TRUE", "FALSE"))
                    )
  )
  bci_onpp <- rbind(bci_onpp, 
                    data.frame(sample = NA,
                               hattheta = NA,
                               lower = NA,
                               upper = NA,
                               include = factor("FALSE", levels = c("TRUE", "FALSE"))
                    )
  )
  bci_onppseq <- rbind(bci_onppseq, 
                       data.frame(sample = NA,
                                  hattheta = NA,
                                  lower = NA,
                                  upper = NA,
                                  include = factor("FALSE", levels = c("TRUE", "FALSE"))
                       )
  )
  
  
  plot_npp <- ggplot() + 
    geom_pointrange(data = bci_npp,
                    mapping = aes(x = sample, y = hattheta,
                                  ymin = lower, ymax = upper, colour = include),
                    na.rm = TRUE,
                    fatten = 0.5) +
    geom_hline(yintercept = true_val, linetype = "longdash") +
    ggtitle(paste(
      "NPP ", "(",
      formatC(coverage$coverage_npp, 
              format = "f", 
              digits = 2),
      ")",
      sep = ""
    )) + 
    theme_gray() +
    scale_color_manual(values = c("TRUE" = "#00BFC4", "FALSE" = "tomato"), drop = FALSE) +
    scale_y_col1 +
    xlab("") +
    ylab("")

  plot_nppseq <- ggplot() +
    geom_pointrange(data = bci_nppseq,
                    mapping = aes(x = sample, y = hattheta,
                                  ymin = lower, ymax = upper, colour = include),
                    na.rm = TRUE,
                    fatten = 0.5) +
    geom_hline(yintercept = true_val, linetype = "longdash") +
    ggtitle(paste(
      "NPP-SEQ ", "(",
      formatC(coverage$coverage_nppseq, 
              format = "f", 
              digits = 2), 
      ")",
      sep = ""
    )) +
    theme_gray() +
    scale_color_manual(values = c("TRUE" = "#00BFC4", "FALSE" = "tomato"), drop = FALSE) +
    scale_y_col2 +
    xlab("") +
    ylab("")

  plot_onpp <- ggplot() +
    geom_pointrange(data = bci_onpp,
                    mapping = aes(x = sample, y = hattheta,
                                  ymin = lower, ymax = upper, colour = include),
                    na.rm = TRUE,
                    fatten = 0.5) +
    geom_hline(yintercept = true_val, linetype = "longdash") +
    ggtitle(paste(
      "ONPP ", "(",
      formatC(coverage$coverage_onpp, 
              format = "f", 
              digits = 2),
      ")",
      sep = ""
    )) +
    theme_gray() +
    scale_color_manual(values = c("TRUE" = "#00BFC4", "FALSE" = "tomato"), drop = FALSE) +
    scale_y_col2 +
    xlab("") +
    ylab("")

  plot_onppseq <- ggplot() +
    geom_pointrange(data = bci_onppseq,
                    mapping = aes(x = sample, y = hattheta,
                                  ymin = lower, ymax = upper, colour = include),
                    na.rm = TRUE,
                    fatten = 0.5) +
    geom_hline(yintercept = true_val, linetype = "longdash") +
    ggtitle(paste(
      "ONPP-SEQ ", "(",
      formatC(coverage$coverage_onppseq, 
              format = "f", 
              digits = 2),
      ")",
      sep = ""
    )) +
    theme_gray() +
    scale_color_manual(values = c("TRUE" = "#00BFC4", "FALSE" = "tomato"), drop = FALSE) +
    scale_y_col2 +
    xlab("") +
    ylab("")

  # combine plots
  title <- paste(prob*100, "%",
                 " BCI for scenario ", 
                 scenario, 
                 " (", 
                 # model, 
                 ")", 
                 sep = "")
  comb_plot <- (plot_npp + plot_nppseq + plot_onpp + plot_onppseq) + 
    plot_layout(guides = "collect", axis_titles = "collect", ncol = 4) &
    theme(legend.position = "none") &
    theme(text = element_text(size = 16),        
          axis.title = element_text(size = 18),  
          axis.text = element_text(size = 16),   
          legend.title = element_text(size = 18),
          legend.text = element_text(size = 16))
  return(comb_plot)
}

plot_sce_lm <- function(j, hat,
                         betastar = get("betastar", envir = .GlobalEnv),
                         sgstar = get("sgstar", envir = .GlobalEnv)) {
  # betastar/sgstar default to looking up globals of the same name (see the
  # note in plot_sce_bin() above) so existing call sites keep working.
  sim <- hat[[j]]
  true_value <- c(betastar, sgstar)
  sce <- as.roman(ceiling(j/3))
  sce <- ifelse(j %% 3 == 1, paste0(sce,".I"), ifelse(j %% 3 == 2, paste0(sce,".II"), paste0(sce,".III")))
  plots_delta <- lapply(seq_along(sim$hatdelta), function(i) {
    plot <- plot_boxplot(sim$hatdelta[[i]], bquote(hat(delta)[.(i)])) + coord_flip() +
      ggtitle(paste("Scenario ", sce))
    if (i>1) {
      plot <- plot + scale_x_discrete(labels = NULL) + ggtitle(NULL)
    }
    return(plot)
  })
  plot_theta <- lapply(seq_along(sim$hattheta), function(i) {
    if(i <= length(betastar)) {
      name <- bquote(hat(beta)[.(i)])
      plot <- plot_boxplot(sim$hattheta[[i]], name) +
        geom_hline(yintercept = true_value[i], linetype = "dotted", color = "red", size = 1) +
        coord_flip()
    } else {
      name <- bquote(hat(sigma^2))
      plot <- plot_boxplot(sim$hattheta[[i]], name) +
        geom_hline(yintercept = true_value[i]^2, linetype = "dotted", color = "red", size = 1) +
        coord_flip()
    }

    if (i>1) {
      plot <- plot + scale_x_discrete(labels = NULL)
    }
    return(plot)
  })
  plot_delta_ <- Reduce(`+`, plots_delta) + plot_layout(ncol = length(plots_delta)) +
    scale_x_discrete(labels = NULL)
  plot_theta_ <- Reduce(`+`, plot_theta) + plot_layout(ncol = length(plot_theta)) +
    scale_x_discrete(labels = NULL)
  return(list(plot_delta_,
              plot_theta_))
}


plot_sim_delta <- function(sce, sim, sim_sces){
  name_sce <- as.roman(ceiling(sce/3))
  name_sce <- ifelse(sce %% 3 == 1, paste0(name_sce,".I"), 
                     ifelse(sce %% 3 == 2, paste0(name_sce,".II"), paste0(name_sce,".III")))
  # plot_delta1 <- ggplot() +
  #   geom_density(aes(x = sim_sces[[sce]]$delta[[sim]]$delta_npp[,1], color = "NPP"), linewidth = 0.8) +
  #   geom_density(aes(x = sim_sces[[sce]]$delta[[sim]]$delta_nppseq[,1], color = "NPP_SEQ"), linewidth = 0.8) +
  #   geom_density(aes(x = sim_sces[[sce]]$delta[[sim]]$delta_onpp[,1], color = "ONPP"), linewidth = 0.8) +
  #   geom_density(aes(x = sim_sces[[sce]]$delta[[sim]]$delta_onppseq[,1], color = "ONPP_SEQ"), linewidth = 0.8) +
  #   scale_color_manual(name = NULL, 
  #                      values = c("NPP" = RColorBrewer::brewer.pal(4, "Set1")[1], 
  #                                 "NPP_SEQ" = RColorBrewer::brewer.pal(4, "Set1")[2],
  #                                 "ONPP" = RColorBrewer::brewer.pal(4, "Set1")[3],
  #                                 "ONPP_SEQ" = RColorBrewer::brewer.pal(4, "Set1")[4])) +
  #   labs(x = expression(delta[1]),
  #        y = NULL) +
  #   ggtitle(paste0("Scenario ", name_sce))
  # plot_delta2 <- ggplot() +
  #   geom_density(aes(x = sim_sces[[sce]]$delta[[sim]]$delta_npp[,2], color = "NPP"), linewidth = 0.8) +
  #   geom_density(aes(x = sim_sces[[sce]]$delta[[sim]]$delta_nppseq[,2], color = "NPP_SEQ"), linewidth = 0.8) +
  #   geom_density(aes(x = sim_sces[[sce]]$delta[[sim]]$delta_onpp[,2], color = "ONPP"), linewidth = 0.8) +
  #   geom_density(aes(x = sim_sces[[sce]]$delta[[sim]]$delta_onppseq[,2], color = "ONPP_SEQ"), linewidth = 0.8) +
  #   scale_color_manual(name = NULL, 
  #                      values = c("NPP" = RColorBrewer::brewer.pal(4, "Set1")[1], 
  #                                 "NPP_SEQ" = RColorBrewer::brewer.pal(4, "Set1")[2],
  #                                 "ONPP" = RColorBrewer::brewer.pal(4, "Set1")[3],
  #                                 "ONPP_SEQ" = RColorBrewer::brewer.pal(4, "Set1")[4])) +
  #   labs(x = expression(delta[2]),
  #        y = NULL)
  # plot_delta3 <- ggplot() +
  #   geom_density(aes(x = sim_sces[[sce]]$delta[[sim]]$delta_npp[,3], color = "NPP"), linewidth = 0.8) +
  #   geom_density(aes(x = sim_sces[[sce]]$delta[[sim]]$delta_nppseq[,3], color = "NPP_SEQ"), linewidth = 0.8) +
  #   geom_density(aes(x = sim_sces[[sce]]$delta[[sim]]$delta_onpp[,3], color = "ONPP"), linewidth = 0.8) +
  #   geom_density(aes(x = sim_sces[[sce]]$delta[[sim]]$delta_onppseq[,3], color = "ONPP_SEQ"), linewidth = 0.8) +
  #   scale_color_manual(name = NULL, 
  #                      values = c("NPP" = RColorBrewer::brewer.pal(4, "Set1")[1], 
  #                                 "NPP_SEQ" = RColorBrewer::brewer.pal(4, "Set1")[2],
  #                                 "ONPP" = RColorBrewer::brewer.pal(4, "Set1")[3],
  #                                 "ONPP_SEQ" = RColorBrewer::brewer.pal(4, "Set1")[4])) +
  #   labs(x = expression(delta[3]),
  #        y = NULL)
  # 
  # combined_plots <- (plot_delta1 + plot_delta2 + plot_delta3 + plot_theta) + 
  #   plot_layout(ncol = 4, guides = "collect")
  
  
  plot_npp <- ggplot() +
    geom_density(aes(x = sim_sces[[sce]]$delta[[sim]]$delta_npp[,1], color = "delta1"), linewidth = 0.8) +
    geom_density(aes(x = sim_sces[[sce]]$delta[[sim]]$delta_npp[,2], color = "delta2"), linewidth = 0.8) +
    geom_density(aes(x = sim_sces[[sce]]$delta[[sim]]$delta_npp[,3], color = "delta3"), linewidth = 0.8) +
    scale_color_manual(name = NULL, 
                       values = c("delta1" = RColorBrewer::brewer.pal(3, "Set1")[1], 
                                  "delta2" = RColorBrewer::brewer.pal(3, "Set1")[2],
                                  "delta3" = RColorBrewer::brewer.pal(3, "Set1")[3]),
                       labels = c(expression(delta[1]), 
                                  expression(delta[2]),
                                  expression(delta[3]))) +
    labs(x = expression(delta),
         y = NULL) +
    ggtitle(label = paste0("Scenario ", name_sce),
            subtitle = "NPP")
  
  plot_nppseq <- ggplot() +
    geom_density(aes(x = sim_sces[[sce]]$delta[[sim]]$delta_nppseq[,1], color = "delta1"), linewidth = 0.8) +
    geom_density(aes(x = sim_sces[[sce]]$delta[[sim]]$delta_nppseq[,2], color = "delta2"), linewidth = 0.8) +
    geom_density(aes(x = sim_sces[[sce]]$delta[[sim]]$delta_nppseq[,3], color = "delta3"), linewidth = 0.8) +
    scale_color_manual(name = NULL, 
                       values = c("delta1" = RColorBrewer::brewer.pal(3, "Set1")[1], 
                                  "delta2" = RColorBrewer::brewer.pal(3, "Set1")[2],
                                  "delta3" = RColorBrewer::brewer.pal(3, "Set1")[3]),
                       labels = c(expression(delta[1]), 
                                  expression(delta[2]),
                                  expression(delta[3]))) +
    labs(x = expression(delta),
         y = NULL) +
    ggtitle(label = NULL,
            subtitle =  "NPP-SEQ")
  
  plot_onpp <- ggplot() +
    geom_density(aes(x = sim_sces[[sce]]$delta[[sim]]$delta_onpp[,1], color = "delta1"), linewidth = 0.8) +
    geom_density(aes(x = sim_sces[[sce]]$delta[[sim]]$delta_onpp[,2], color = "delta2"), linewidth = 0.8) +
    geom_density(aes(x = sim_sces[[sce]]$delta[[sim]]$delta_onpp[,3], color = "delta3"), linewidth = 0.8) +
    scale_color_manual(name = NULL, 
                       values = c("delta1" = RColorBrewer::brewer.pal(3, "Set1")[1], 
                                  "delta2" = RColorBrewer::brewer.pal(3, "Set1")[2],
                                  "delta3" = RColorBrewer::brewer.pal(3, "Set1")[3]),
                       labels = c(expression(delta[1]), 
                                  expression(delta[2]),
                                  expression(delta[3]))) +
    labs(x = expression(delta),
         y = NULL) +
    ggtitle(label = NULL,
            subtitle = "ONPP")
  
  plot_onppseq <- ggplot() +
    geom_density(aes(x = sim_sces[[sce]]$delta[[sim]]$delta_onppseq[,1], color = "delta1"), linewidth = 0.8) +
    geom_density(aes(x = sim_sces[[sce]]$delta[[sim]]$delta_onppseq[,2], color = "delta2"), linewidth = 0.8) +
    geom_density(aes(x = sim_sces[[sce]]$delta[[sim]]$delta_onppseq[,3], color = "delta3"), linewidth = 0.8) +
    scale_color_manual(name = NULL, 
                       values = c("delta1" = RColorBrewer::brewer.pal(3, "Set1")[1], 
                                  "delta2" = RColorBrewer::brewer.pal(3, "Set1")[2],
                                  "delta3" = RColorBrewer::brewer.pal(3, "Set1")[3]),
                       labels = c(expression(delta[1]), 
                                  expression(delta[2]),
                                  expression(delta[3]))) +
    labs(x = expression(delta),
         y = NULL) +
    ggtitle(label = NULL,
            subtitle = "ONPP-SEQ")
  combined_plots <- (plot_npp + plot_nppseq + plot_onpp + plot_onppseq) + 
    plot_layout(ncol = 4, guides = "collect")
  return(combined_plots)
}

plot_sim_theta <- function(sce, sim, sim_sces, true_value = get("true_value", envir = .GlobalEnv)) {
  # true_value defaults to looking up the global of the same name (see the
  # note in plot_sce_bin() above) so existing call sites keep working.
  name_sce <- as.roman(ceiling(sce/3))
  name_sce <- ifelse(sce %% 3 == 1, paste0(name_sce,".I"),
                     ifelse(sce %% 3 == 2, paste0(name_sce,".II"), paste0(name_sce,".III")))
  plot_theta <- ggplot() +
    geom_density(aes(x = sim_sces[[sce]]$theta[[sim]]$theta_npp, color = "NPP"), linewidth = 0.8) +
    geom_density(aes(x = sim_sces[[sce]]$theta[[sim]]$theta_nppseq, color = "NPP_SEQ"), linewidth = 0.8) +
    geom_density(aes(x = sim_sces[[sce]]$theta[[sim]]$theta_onpp, color = "ONPP"), linewidth = 0.8) +
    geom_density(aes(x = sim_sces[[sce]]$theta[[sim]]$theta_onppseq, color = "ONPP_SEQ"), linewidth = 0.8) +
    geom_vline(xintercept = true_value, linetype = "dotted", color = "black", size = 1) +
    scale_color_manual(name = NULL,
                       values = c("NPP" = RColorBrewer::brewer.pal(4, "Set1")[1],
                                  "NPP_SEQ" = RColorBrewer::brewer.pal(4, "Set1")[2],
                                  "ONPP" = RColorBrewer::brewer.pal(4, "Set1")[3],
                                  "ONPP_SEQ" = RColorBrewer::brewer.pal(4, "Set1")[4])) +
    labs(x = expression(theta),
         y = NULL) +
    xlim(0,1) +
    ggtitle(label = paste0("Scenario ", name_sce))
  return(plot_theta)
}

# ------------------------------------------------------------------------------
# Shared look for the sce3-across-priors figures
# ------------------------------------------------------------------------------
# Both plot_sce3_theta()/plot_sce3_eta() used to build one ggplot per panel
# by hand (one lapply per row), set a title/subtitle on each individually,
# and stitch the panels together with patchwork's wrap_plots(). That meant
# every row had to remember, on its own, to add the "(A)"/"(B)"/"(C)"/"(D)"
# column labels -- plot_sce3_eta()'s NPP row did, but its ONPP row never
# did, so reading an ONPP panel meant counting columns against the row
# above it to figure out which prior it belonged to. Both functions are
# rewritten below to reshape their data into one long data.frame and use
# facet_grid()/facet_wrap() instead -- ggplot then generates consistent
# strip labels for every row for free, and there's only one plot object
# (and one legend) per scenario instead of 4-8 hand-assembled ones.
#
# Okabe-Ito-derived, colorblind-safe palette (avoids the red/green pairing
# in the old RColorBrewer "Set1" palette, which is exactly the confusable
# pair for the most common form of color blindness) -- linetype is mapped
# redundantly alongside color so the lines stay distinguishable in
# grayscale print, too.
.sce3_theme <- function() {
  theme_bw(base_size = 11) +
    theme(
      strip.background  = element_rect(fill = "grey90", color = NA),
      strip.text         = element_text(face = "bold"),
      plot.title         = element_text(face = "bold"),
      panel.grid.minor   = element_blank(),
      legend.position    = "bottom"
    )
}

plot_sce3_theta <- function(sub_sce, sim, sim_sces, true_value = get("true_value", envir = .GlobalEnv),
                            prior_labels = c("(A) Default", "(B) \u03B7 flat",
                                             "(C) \u03B7, \u03B3 flat", "(D) All flat"),
                            # The 3 sub-scenarios of Scenario III (incompatible order) in the
                            # bernoulli study, in order -- same theta0/Discrepancy values as
                            # analysis_bin.r's df_mse_theta rows 7-9 (scenarios_bin.R's
                            # sce3_theta0). Order is fixed ("Incompatible") since plot_sce3_*()
                            # only ever plots Scenario III's sub-scenarios. No "Scenario III.I"
                            # label -- theta0/Order/Discrepancy already identify the scenario,
                            # matching the analysis_seq_*.r figures' title pattern.
                            sce_order = "Incompatible",
                            sce_discrepancy = c("None", "Small", "Large"),
                            sce_theta0 = list(c(0.50, 0.40, 0.30), c(0.40, 0.30, 0.20), c(0.25, 0.15, 0.05))) {
  # true_value defaults to looking up the global of the same name (see the
  # note in plot_sce_bin() above) so existing call sites keep working.
  theta0_str <- paste0("(", paste(formatC(sce_theta0[[sub_sce]], format = "f", digits = 2), collapse = ", "), ")")
  full_title <- bquote(
    theta[0] == .(theta0_str) * ", Order: " * .(sce_order) * ", Discrepancy: " * .(sce_discrepancy[sub_sce])
  )
  sim_sub_sce <- sim_sces[[sub_sce]]

  # th$theta_npp/theta_onpp are `posterior` draws objects, not plain
  # numeric vectors -- rbind-ing data.frame()s built directly from them
  # trips base R's match.names() (draws-object metadata riding along on
  # each column), the same failure already fixed the same way in
  # plot_eta_theta_density(). as.numeric() strips that metadata first.
  long_df <- do.call(rbind, lapply(seq_along(sim_sub_sce), function(i) {
    th <- sim_sub_sce[[i]]$theta[[sim]]
    rbind(
      data.frame(value = as.numeric(th$theta_npp),  method = "NPP",  prior = prior_labels[i]),
      data.frame(value = as.numeric(th$theta_onpp), method = "ONPP", prior = prior_labels[i])
    )
  }))
  long_df$method <- factor(long_df$method, levels = c("NPP", "ONPP"))
  long_df$prior  <- factor(long_df$prior, levels = prior_labels)

  ggplot(long_df, aes(x = value, color = method)) +
    geom_density(linewidth = 0.8) +
    geom_vline(xintercept = true_value, linetype = "dotted", color = "black", linewidth = 0.6) +
    facet_wrap(~ prior, nrow = 1, scales = "free_y") +
    # Okabe-Ito blue/vermillion -- deliberately NOT reusing any of the eta
    # plot's 3 colors (orange/sky-blue/bluish-green, see plot_sce3_eta()
    # below), since these are companion figures for the same scenario and
    # the same color meaning two different things across them is confusing.
    scale_color_manual(name = NULL, values = c("NPP" = "#0072B2", "ONPP" = "#D55E00")) +
    scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
    labs(x = expression(theta), y = NULL, title = full_title) +
    .sce3_theme()
}

plot_sce3_eta <- function(sub_sce, sim, sim_sces,
                          prior_labels = c("(A) Default", "(B) \u03B7 flat",
                                           "(C) \u03B7, \u03B3 flat", "(D) All flat"),
                          # Same rationale as plot_sce3_theta() above.
                          sce_order = "Incompatible",
                          sce_discrepancy = c("None", "Small", "Large"),
                          sce_theta0 = list(c(0.50, 0.40, 0.30), c(0.40, 0.30, 0.20), c(0.25, 0.15, 0.05))) {
  theta0_str <- paste0("(", paste(formatC(sce_theta0[[sub_sce]], format = "f", digits = 2), collapse = ", "), ")")
  full_title <- bquote(
    theta[0] == .(theta0_str) * ", Order: " * .(sce_order) * ", Discrepancy: " * .(sce_discrepancy[sub_sce])
  )
  sim_sub_sce <- sim_sces[[sub_sce]]

  # Same fix as plot_sce3_theta() above: delta_npp/delta_onpp columns carry
  # `posterior` draws-object metadata, which breaks rbind()'s match.names()
  # once several of these data.frames get stacked -- as.numeric() strips it.
  long_df <- do.call(rbind, lapply(seq_along(sim_sub_sce), function(i) {
    delta <- sim_sub_sce[[i]]$delta[[sim]]
    rbind(
      data.frame(value = as.numeric(delta$delta_npp[, 1]),  component = "eta[1]", method = "NPP",  prior = prior_labels[i]),
      data.frame(value = as.numeric(delta$delta_npp[, 2]),  component = "eta[2]", method = "NPP",  prior = prior_labels[i]),
      data.frame(value = as.numeric(delta$delta_npp[, 3]),  component = "eta[3]", method = "NPP",  prior = prior_labels[i]),
      data.frame(value = as.numeric(delta$delta_onpp[, 1]), component = "eta[1]", method = "ONPP", prior = prior_labels[i]),
      data.frame(value = as.numeric(delta$delta_onpp[, 2]), component = "eta[2]", method = "ONPP", prior = prior_labels[i]),
      data.frame(value = as.numeric(delta$delta_onpp[, 3]), component = "eta[3]", method = "ONPP", prior = prior_labels[i])
    )
  }))
  long_df$method    <- factor(long_df$method, levels = c("NPP", "ONPP"))
  long_df$prior     <- factor(long_df$prior, levels = prior_labels)
  long_df$component <- factor(long_df$component, levels = c("eta[1]", "eta[2]", "eta[3]"))

  ggplot(long_df, aes(x = value, color = component)) +
    geom_density(bounds = c(0, 1), linewidth = 0.8) +
    # Plain facet_grid(scales = "free_y") only frees the y-axis per ROW
    # (shared across all 4 priors in that row), and facet_wrap() frees it
    # per panel but loses the shared prior/method header (each panel gets
    # its own repeated "NPP"/"(A) Default" strip instead). ggh4x's
    # facet_grid2(independent = "y") keeps facet_grid's shared header
    # (priors across the top, method down the side, each once) while still
    # giving every individual panel its own free y-axis -- requires the
    # ggh4x package (library(ggh4x) at the top of the calling script).
    ggh4x::facet_grid2(method ~ prior, scales = "free_y", independent = "y") +
    scale_color_manual(
      name   = NULL,
      breaks = c("eta[1]", "eta[2]", "eta[3]"),
      labels = c(expression(eta[1]), expression(eta[2]), expression(eta[3])),
      values = c("eta[1]" = "#E69F00", "eta[2]" = "#56B4E9", "eta[3]" = "#009E73")
    ) +
    scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
    labs(x = expression(eta), y = NULL, title = full_title) +
    .sce3_theme()
}

# ------------------------------------------------------------------------------
# Shared Default-vs-SEQ prior/posterior density plots
# ------------------------------------------------------------------------------
#
# analysis_seq_bin.r, analysis_seq_normal.r and analysis_seq_lm.r each used
# to carry their own copy of plot_eta_bin()/plot_eta_normal()/plot_eta_lm()
# -- three byte-for-byte identical implementations (same delta_npp/
# delta_nppseq/delta_onpp/delta_onppseq structure, same 3 eta components,
# same 0-1 domain) built out of 3 hand-assembled ggplot() calls stitched
# together with patchwork, each repeating its own scale_color_manual() and
# theme_minimal(). That's exactly the copy-paste-across-families pattern
# this file's other comments call out as the source of the LM
# hist_1/hist_3 relabeling bug -- one shared function here instead, reshaped
# to a long data.frame and faceted the same way as plot_sce3_eta() above,
# so there's one definition to fix instead of three that can drift apart.
#
# Default/SEQ get their own pair -- black and Okabe-Ito reddish-purple --
# instead of reusing plot_sce3_theta()'s NPP/ONPP blue/vermillion: the two
# figures compare different things (Default-vs-SEQ here, NPP-vs-ONPP
# there), so the same colors meaning something different across them is
# confusing, exactly like the theta/eta color clash fixed earlier in this
# file. Defined once here so plot_eta_prior_post(), plot_theta_prior_post()
# and analysis_seq_lm.r's plot_theta_lm() (whose theta is a 4-vector and
# stays separate, but shares this same Default/SEQ legend) can't drift
# apart the same way the three old per-family copies did.
.default_seq_colors <- c("Default" = "#000000", "SEQ" = "#CC79A7")

plot_eta_prior_post <- function(draws, ordered = FALSE, post = FALSE, title = NULL, leg = TRUE) {
  eta_nppseq <- draws$delta_nppseq
  eta_npp    <- draws$delta_npp
  if (ordered) {
    eta_nppseq <- draws$delta_onppseq
    eta_npp    <- draws$delta_onpp
  }

  # as.numeric(): see the match.names() note on plot_sce3_eta() above --
  # same posterior/draws-object metadata issue when several of these get
  # combined via rbind().
  long_df <- do.call(rbind, lapply(1:3, function(k) {
    rbind(
      data.frame(value = as.numeric(eta_npp[, k]),    method = "Default", component = paste0("eta[", k, "]")),
      data.frame(value = as.numeric(eta_nppseq[, k]), method = "SEQ",     component = paste0("eta[", k, "]"))
    )
  }))
  long_df$method    <- factor(long_df$method, levels = c("Default", "SEQ"))
  long_df$component <- factor(long_df$component, levels = paste0("eta[", 1:3, "]"))

  y_lab <- if (post) {
    expression(paste(pi[A], "(", eta ~ "|" ~ D ~ "," ~ D[0], ")"))
  } else {
    expression(paste(pi[A], "(", eta, ")"))
  }

  p <- ggplot(long_df, aes(x = value, color = method)) +
    # Prior curves dotted, posterior curves solid -- same linetype, on top
    # of the existing color coding, so a prior/posterior pair is
    # distinguishable even if printed in grayscale or the colors are hard
    # to tell apart.
    geom_density(bounds = c(0, 1), linewidth = 0.8, linetype = if (post) "solid" else "dashed") +
    facet_wrap(~ component, nrow = 1, labeller = label_parsed, scales = "free_y") +
    scale_color_manual(name = NULL, values = .default_seq_colors) +
    scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
    labs(x = expression(eta), y = y_lab, title = title) +
    .sce3_theme()

  p + theme(legend.position = if (leg) "bottom" else "none")
}

# Shared scalar-theta counterpart, for the two families (bin, normal) whose
# theta is a single number -- lm's theta is a 4-vector (beta0/beta1/beta2/
# sigma^2) with its own MLE-marker overlay, different enough that
# plot_theta_lm() stays a separate, hand-written function in
# analysis_seq_lm.r rather than being forced into this shape.
plot_theta_prior_post <- function(draws, ordered = FALSE, post = FALSE, theta_true, low_band, upp_band,
                                  title = NULL, leg = TRUE) {
  theta_nppseq <- draws$theta_nppseq
  theta_npp    <- draws$theta_npp
  if (ordered) {
    theta_nppseq <- draws$theta_onppseq
    theta_npp    <- draws$theta_onpp
  }

  long_df <- rbind(
    data.frame(value = as.numeric(theta_npp),    method = "Default"),
    data.frame(value = as.numeric(theta_nppseq), method = "SEQ")
  )
  long_df$method <- factor(long_df$method, levels = c("Default", "SEQ"))

  y_lab <- if (post) {
    expression(paste(pi[0], "(", theta ~ "|" ~ D ~ "," ~ D[0], ")"))
  } else {
    expression(paste(pi[0], "(", theta, ")"))
  }

  p <- ggplot(long_df, aes(x = value, color = method)) +
    # Prior curves dotted, posterior curves solid -- same reasoning as
    # plot_eta_prior_post() above.
    geom_density(linewidth = 0.8, bounds = c(theta_true - low_band, theta_true + upp_band),
                linetype = if (post) "solid" else "dashed") +
    # Reference line uses grey, not black, so it doesn't visually collide
    # with the (now black) "Default" density curve.
    geom_vline(xintercept = theta_true, linetype = "dotted", color = "grey40", linewidth = 0.6) +
    scale_color_manual(name = NULL, values = .default_seq_colors) +
    scale_x_continuous(limits = c(theta_true - low_band, theta_true + upp_band)) +
    labs(x = expression(theta), y = y_lab, title = title) +
    .sce3_theme()

  p + theme(legend.position = if (leg) "bottom" else "none")
}

plot_eta_prior <- function(samples){
  plots_npp <- lapply(1:length(samples), function(i) {
    plot <- ggplot() +
      geom_density(bounds = c(0, 1), aes(x = samples[[i]]$delta_npp[,1], color = "delta1"), linewidth = 0.8) +
      geom_density(bounds = c(0, 1), aes(x = samples[[i]]$delta_npp[,2], color = "delta2"), linewidth = 0.8) +
      geom_density(bounds = c(0, 1), aes(x = samples[[i]]$delta_npp[,3], color = "delta3"), linewidth = 0.8) +
      scale_color_manual(name = NULL, 
                         values = c("delta1" = RColorBrewer::brewer.pal(3, "Set1")[1], 
                                    "delta2" = RColorBrewer::brewer.pal(3, "Set1")[2],
                                    "delta3" = RColorBrewer::brewer.pal(3, "Set1")[3]),
                         labels = c(expression(eta[1]), 
                                    expression(eta[2]),
                                    expression(eta[3]))) +
      labs(x = expression(eta),
           y = NULL)
    if (i == 1) {
      plot <- plot +
          ggtitle(label = "NPP",
                  subtitle = paste0("(", LETTERS[i], ")"))
    } else {
      plot <- plot + 
        ggtitle(label = NULL,
                subtitle = paste0("(", LETTERS[i], ")"))
    }
    return(plot)
    })
  
  plots_onpp <- lapply(1:length(samples), function(i) {
    plot <- ggplot() +
      geom_density(bounds = c(0, 1), aes(x = samples[[i]]$delta_onpp[,1], color = "delta1"), linewidth = 0.8) +
      geom_density(bounds = c(0, 1), aes(x = samples[[i]]$delta_onpp[,2], color = "delta2"), linewidth = 0.8) +
      geom_density(bounds = c(0, 1), aes(x = samples[[i]]$delta_onpp[,3], color = "delta3"), linewidth = 0.8) +
      scale_color_manual(name = NULL, 
                         values = c("delta1" = RColorBrewer::brewer.pal(3, "Set1")[1], 
                                    "delta2" = RColorBrewer::brewer.pal(3, "Set1")[2],
                                    "delta3" = RColorBrewer::brewer.pal(3, "Set1")[3]),
                         labels = c(expression(eta[1]), 
                                    expression(eta[2]),
                                    expression(eta[3]))) +
      labs(x = expression(eta),
           y = NULL)
      if (i == 1) {
        plot <- plot +
          ggtitle(label = "ONPP")
      }
    return(plot)
  })

  plots_onpp_gamma <- lapply(1:length(samples), function(i) {
    plot <- ggplot() +
      geom_density(bounds = c(0, 1), aes(x = samples[[i]]$delta_onpp[,1], color = "delta1"), linewidth = 0.8) +
      geom_density(bounds = c(0, 1), aes(x = samples[[i]]$delta_onpp[,2] - samples[[i]]$delta_onpp[,1], 
        color = "delta2"), linewidth = 0.8) +
      geom_density(aes(x = samples[[i]]$delta_onpp[,3] - samples[[i]]$delta_onpp[,2], 
        color = "delta3"), linewidth = 0.8) +
      scale_color_manual(name = NULL, 
                         values = c("delta1" = RColorBrewer::brewer.pal(3, "Set3")[1], 
                                    "delta2" = RColorBrewer::brewer.pal(3, "Set3")[3],
                                    "delta3" = RColorBrewer::brewer.pal(3, "Set3")[5]),
                         labels = c(expression(gamma[1]), 
                                    expression(gamma[2]),
                                    expression(gamma[3]))) +
      labs(x = expression(gamma),
           y = NULL)
    return(plot)
  })
  
  combined_plots <- wrap_plots(plots_npp, ncol = length(plots_npp)) / 
    wrap_plots(plots_onpp, ncol = length(plots_onpp)) /
    wrap_plots(plots_onpp_gamma, ncol = length(plots_onpp_gamma)) +
    plot_layout(guides = "collect")
  return(combined_plots)
}

plot_theta_prior <- function(samples){
  plot_theta <- lapply(1:length(samples), function(i) {
    plot <- 
    ggplot() +
    geom_density(aes(x = samples[[i]], color = "dens"), linewidth = 0.8) +
    scale_color_manual(name = NULL,
                       values = c("dens" = RColorBrewer::brewer.pal(4, "Set2")[3])) +
    labs(x = expression(theta),
         y = NULL) +
    xlim(0,1) + 
      theme(legend.position = "none") +
        ggtitle(label = NULL,
                subtitle = paste0("(", LETTERS[i], ")"))
    return(plot)
  })
  combined_plots <- wrap_plots(plot_theta, ncol = length(plot_theta)) + 
    plot_layout(guides = "collect")
  return(combined_plots)
  }


plot_eta_post <- function(samples, sim, seq_norm = F){
  title_npp <- paste0("NPP")
  title_onpp <- paste0("ONPP")
  if (seq_norm) {
    samples$delta[[sim]]$delta_npp <- samples$delta[[sim]]$delta_nppseq
    samples$delta[[sim]]$delta_onpp <- samples$delta[[sim]]$delta_onppseq
    title_npp <- paste0("NPP-SEQ")
    title_onpp <- paste0("ONPP-SEQ")
  }
  plot_npp <- ggplot() +
    geom_density(bounds = c(0, 1), aes(x = samples$delta[[sim]]$delta_npp[,1], color = "delta1"), linewidth = 0.8) +
    geom_density(bounds = c(0, 1), aes(x = samples$delta[[sim]]$delta_npp[,2], color = "delta2"), linewidth = 0.8) +
    geom_density(bounds = c(0, 1), aes(x = samples$delta[[sim]]$delta_npp[,3], color = "delta3"), linewidth = 0.8) +
    scale_color_manual(name = NULL, 
                        values = c("delta1" = RColorBrewer::brewer.pal(3, "Set1")[1], 
                                  "delta2" = RColorBrewer::brewer.pal(3, "Set1")[2],
                                  "delta3" = RColorBrewer::brewer.pal(3, "Set1")[3]),
                        labels = c(expression(eta[1]), 
                                  expression(eta[2]),
                                  expression(eta[3]))) +
    labs(x = expression(eta),
          y = NULL) +
    ggtitle(label = title_npp)
  
  plot_onpp <- ggplot() +
    geom_density(bounds = c(0, 1), aes(x = samples$delta[[sim]]$delta_onpp[,1], color = "delta1"), linewidth = 0.8) +
    geom_density(bounds = c(0, 1), aes(x = samples$delta[[sim]]$delta_onpp[,2], color = "delta2"), linewidth = 0.8) +
    geom_density(bounds = c(0, 1), aes(x = samples$delta[[sim]]$delta_onpp[,3], color = "delta3"), linewidth = 0.8) +
    scale_color_manual(name = NULL, 
                        values = c("delta1" = RColorBrewer::brewer.pal(3, "Set1")[1], 
                                  "delta2" = RColorBrewer::brewer.pal(3, "Set1")[2],
                                  "delta3" = RColorBrewer::brewer.pal(3, "Set1")[3]),
                        labels = c(expression(eta[1]), 
                                  expression(eta[2]),
                                  expression(eta[3]))) +
    labs(x = expression(eta),
          y = NULL) +
    ggtitle(label = title_onpp)

  plot_onpp_gamma <- ggplot() +
    geom_density(bounds = c(0, 1), aes(x = samples$delta[[sim]]$delta_onpp[,1], color = "delta1"), linewidth = 0.8) +
    geom_density(bounds = c(0, 1), aes(x = samples$delta[[sim]]$delta_onpp[,2] - samples$delta[[sim]]$delta_onpp[,1], 
      color = "delta2"), linewidth = 0.8) +
    geom_density(aes(x = samples$delta[[sim]]$delta_onpp[,3] - samples$delta[[sim]]$delta_onpp[,2], 
      color = "delta3"), linewidth = 0.8) +
    scale_color_manual(name = NULL, 
                        values = c("delta1" = RColorBrewer::brewer.pal(3, "Set3")[1], 
                                  "delta2" = RColorBrewer::brewer.pal(3, "Set3")[3],
                                  "delta3" = RColorBrewer::brewer.pal(3, "Set3")[5]),
                        labels = c(expression(gamma[1]), 
                                  expression(gamma[2]),
                                  expression(gamma[3]))) +
    labs(x = expression(gamma),
          y = NULL)
  
  combined_plots <- plot_npp + plot_onpp + plot_onpp_gamma +
    plot_layout(guides = "collect")
  return(combined_plots)
}

# =============================================================================
# General, family-agnostic density plot for a single (scenario, replicate)
# =============================================================================
#
# All four families -- bin, poi, normal (normal_fixed_var), lm -- funnel
# through .assemble_sce_result() and sim_sce() into the exact same shape for
# a single scenario's simulation object `sim`:
#
#   sim$delta[[replicate]]$delta_npp / delta_nppseq / delta_onpp / delta_onppseq
#   sim$theta[[replicate]]$theta_npp / theta_nppseq / theta_onpp / theta_onppseq
#
# so this function works unchanged across families: `family` is used only
# for axis bounds/labels, never to decide how the data is read. The one real
# structural difference -- theta being a single vector of draws (bin/poi/
# normal, where theta is a scalar parameter) vs. a matrix (lm, one column
# per beta + sigma) -- is detected automatically with is.matrix(), not
# hardcoded per family.
#
# Arguments:
#   sim            one scenario's simulation object, e.g. sim_sces[[7]] for
#                   bin/normal, or sim_sces_beta_n0_ge_n[[2]] for lm -- i.e.
#                   whatever list of scenarios that family's analysis script
#                   already loads. Pass the *scenario*, not the whole list.
#   replicate      integer index into sim$delta / sim$theta picking which
#                   simulated replicate to plot (1..num_sim)
#   family         one of "bin", "poi", "normal", "lm" -- only affects
#                   default axis bounds/labels, see above
#   true_theta     optional true theta value(s) -- a single number for
#                   bin/poi/normal, or a vector matching the columns of
#                   theta (e.g. c(betastar, sgstar) for lm)
#   theta_labels   optional list of expressions naming each theta component;
#                   defaults to expression(theta) when theta is univariate,
#                   or beta[1], beta[2], ..., sigma when it's a matrix (lm)
#   scenario_label optional string (e.g. "Scenario III.II") shown in the title
#
# Returns a single patchwork object: one eta density panel per method
# (NPP / NPP-SEQ / ONPP / ONPP-SEQ) on top, theta density panel(s) below --
# for lm, one theta panel per beta/sigma component, each overlaying the four
# methods; for bin/poi/normal, a single panel overlaying the four methods,
# same as plot_sim_theta().
#
# Examples:
#   # bin, scenario III.II (7th element), replicate 3:
#   plot_eta_theta_density(sim_sces[[7]], replicate = 3, family = "bin",
#                           true_theta = 0.5,
#                           scenario_label = "Scenario III.II")
#
#   # lm, "varying beta", n0 >= n, scenario 2, replicate 10:
#   plot_eta_theta_density(sim_sces_beta_n0_ge_n[[2]], replicate = 10,
#                           family = "lm", true_theta = c(betastar, sgstar))
plot_eta_theta_density <- function(sim, replicate,
                                    family = c("bin", "poi", "normal", "lm"),
                                    true_theta = NULL,
                                    theta_labels = NULL, scenario_label = NULL) {
  family <- match.arg(family)
  methods <- c("npp", "nppseq", "onpp", "onppseq")
  method_titles <- c(npp = "NPP", nppseq = "NPP-SEQ", onpp = "ONPP", onppseq = "ONPP-SEQ")

  if (is.null(sim$delta) || is.null(sim$theta)) {
    stop("`sim` doesn't look like a single scenario's simulation object -- ",
         "it has no $delta/$theta. Did you pass the whole sim_sces list ",
         "instead of one element of it (e.g. sim_sces[[sce]])?")
  }
  if (replicate < 1 || replicate > length(sim$delta)) {
    stop("replicate ", replicate, " out of range 1:", length(sim$delta))
  }

  eta_rec <- sim$delta[[replicate]]
  theta_rec <- sim$theta[[replicate]]

  # ---- eta density panels: one per method, all eta components overlaid ----
  n_eta <- ncol(eta_rec$delta_npp)
  eta_pal <- RColorBrewer::brewer.pal(max(3, n_eta), "Set1")[seq_len(n_eta)]
  eta_key <- paste0("eta", seq_len(n_eta))
  eta_expr <- lapply(seq_len(n_eta), function(i) bquote(eta[.(i)]))

  eta_plots <- lapply(methods, function(m) {
    draws <- eta_rec[[paste0("delta_", m)]]
    df <- do.call(rbind, lapply(seq_len(n_eta), function(i) {
      # as.numeric() strips the posterior package's draws_matrix class --
      # draws[, i] otherwise stays a 1-column draws_matrix (its `[` method
      # doesn't drop dimensions), which carries its own as.data.frame()
      # method that names the column after the Stan variable (e.g.
      # "eta[1]") instead of "value". With three components pulled from
      # possibly different draws_matrix objects, that produced three
      # differently-named single-column data.frames and rbind() failed
      # with "names do not match previous names".
      data.frame(value = as.numeric(draws[, i]), component = eta_key[i])
    }))
    p <- ggplot(df) +
      geom_density(bounds = c(0, 1), aes(x = value, color = component), linewidth = 0.8) +
      scale_color_manual(name = NULL,
                         values = setNames(eta_pal, eta_key),
                         labels = eta_expr) +
      labs(x = expression(eta), y = NULL) +
      ggtitle(method_titles[[m]])
    p
  })
  eta_combined <- wrap_plots(eta_plots, ncol = length(eta_plots)) +
    plot_layout(guides = "collect")

  # ---- theta density panel(s) ----------------------------------------------
  is_multivariate_theta <- is.matrix(theta_rec$theta_npp)
  n_theta <- if (is_multivariate_theta) ncol(theta_rec$theta_npp) else 1

  if (is.null(theta_labels)) {
    theta_labels <- if (!is_multivariate_theta) {
      list(expression(theta))
    } else {
      c(lapply(seq_len(n_theta - 1), function(i) bquote(beta[.(i)])),
        list(expression(sigma)))
    }
  }

  method_pal <- setNames(RColorBrewer::brewer.pal(4, "Set1"), unname(method_titles))

  theta_plots <- lapply(seq_len(n_theta), function(j) {
    df <- do.call(rbind, lapply(methods, function(m) {
      draws <- theta_rec[[paste0("theta_", m)]]
      # Same as.numeric() fix as the eta loop above, for the same reason.
      vals <- if (is_multivariate_theta) as.numeric(draws[, j]) else as.numeric(draws)
      data.frame(value = vals, method = method_titles[[m]])
    }))
    p <- ggplot(df) +
      geom_density(aes(x = value, color = method), linewidth = 0.8) +
      scale_color_manual(name = NULL, values = method_pal) +
      labs(x = theta_labels[[j]], y = NULL)
    if (family == "bin") p <- p + xlim(0, 1)  # theta is a probability only for bin
    if (!is.null(true_theta)) {
      true_theta_j <- if (is_multivariate_theta) true_theta[j] else true_theta[1]
      p <- p + geom_vline(xintercept = true_theta_j, linetype = "dotted",
                          color = "black", linewidth = 0.6)
    }
    p
  })
  theta_combined <- wrap_plots(theta_plots, ncol = n_theta) +
    plot_layout(guides = "collect")

  title <- paste0(
    if (!is.null(scenario_label)) paste0(scenario_label, " — ") else "",
    "Replicate ", replicate, " (", toupper(family), ")"
  )

  (eta_combined / theta_combined) +
    plot_annotation(title = title)
}

# =============================================================================
# Same eta/theta density plot as plot_eta_theta_density() above, but in the
# visual style used elsewhere in this file for the sce3/prior-post figures
# (facet by eta component instead of one ggplot per method, .sce3_theme(),
# an Okabe-Ito colorblind-safe palette, a single legend collected at the
# bottom) instead of plot_eta_theta_density()'s one-panel-per-method/
# theme_grey defaults. Family-agnostic like plot_eta_theta_density(); takes
# `title` directly (a string or bquote()/expression, e.g. built the same
# way as analysis_seq_*.r's theta[0] == ...~", Order:"~... titles) instead
# of assembling one from a scenario_label + replicate + family, since
# callers of this one already build that title themselves.
# =============================================================================
plot_eta_theta_density_styled <- function(sim, replicate,
                                           family = c("bin", "poi", "normal", "lm"),
                                           true_theta = NULL,
                                           theta_labels = NULL, title = NULL) {
  family <- match.arg(family)
  methods <- c("npp", "nppseq", "onpp", "onppseq")
  method_labels <- c(npp = "NPP", nppseq = "NPP-SEQ", onpp = "ONPP", onppseq = "ONPP-SEQ")
  # NPP/NPP-SEQ reuse the black/reddish-purple "Default"/"SEQ" pair used
  # throughout this file (see .default_seq_colors); ONPP/ONPP-SEQ get the
  # same treatment in blue/vermillion, so the two order-of-borrowing
  # families read as two paired color groups instead of four arbitrary
  # colors.
  method_colors <- c(NPP = "#000000", "NPP-SEQ" = "#CC79A7",
                     ONPP = "#0072B2", "ONPP-SEQ" = "#D55E00")

  if (is.null(sim$delta) || is.null(sim$theta)) {
    stop("`sim` doesn't look like a single scenario's simulation object -- ",
         "it has no $delta/$theta. Did you pass the whole sim_sces list ",
         "instead of one element of it (e.g. sim_sces[[sce]])?")
  }
  if (replicate < 1 || replicate > length(sim$delta)) {
    stop("replicate ", replicate, " out of range 1:", length(sim$delta))
  }

  eta_rec <- sim$delta[[replicate]]
  theta_rec <- sim$theta[[replicate]]

  # ---- eta: facet by method, color = component -- same component palette
  # (orange/sky-blue/bluish-green) and color-by-component idea as
  # plot_sce3_eta() above, just faceted on method (one scenario here, so
  # there's no "prior" dimension to also facet on like that one does).
  n_eta <- ncol(eta_rec$delta_npp)
  eta_key  <- paste0("eta[", seq_len(n_eta), "]")
  eta_expr <- lapply(seq_len(n_eta), function(i) bquote(eta[.(i)]))
  eta_pal  <- c("eta[1]" = "#E69F00", "eta[2]" = "#56B4E9", "eta[3]" = "#009E73")

  eta_long <- do.call(rbind, lapply(methods, function(m) {
    draws <- eta_rec[[paste0("delta_", m)]]
    do.call(rbind, lapply(seq_len(n_eta), function(i) {
      # as.numeric(): see the match.names() note on plot_sce3_eta() above --
      # same posterior-package draws-object metadata issue.
      data.frame(value = as.numeric(draws[, i]),
                 component = eta_key[i],
                 method = method_labels[[m]])
    }))
  }))
  eta_long$method    <- factor(eta_long$method, levels = unname(method_labels))
  eta_long$component <- factor(eta_long$component, levels = eta_key)

  eta_plot <- ggplot(eta_long, aes(x = value, color = component)) +
    geom_density(bounds = c(0, 1), linewidth = 0.8) +
    facet_wrap(~ method, nrow = 1, scales = "free_y") +
    scale_color_manual(name = NULL, breaks = eta_key, labels = eta_expr, values = eta_pal) +
    scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
    labs(x = expression(eta), y = NULL, title = title) +
    .sce3_theme()

  # ---- theta: single panel (bin/poi/normal) or one per component (lm) ----
  is_multivariate_theta <- is.matrix(theta_rec$theta_npp)
  n_theta <- if (is_multivariate_theta) ncol(theta_rec$theta_npp) else 1

  if (is.null(theta_labels)) {
    theta_labels <- if (!is_multivariate_theta) {
      list(expression(theta))
    } else {
      c(lapply(seq_len(n_theta - 1), function(i) bquote(beta[.(i)])),
        list(expression(sigma)))
    }
  }

  theta_plots <- lapply(seq_len(n_theta), function(j) {
    df <- do.call(rbind, lapply(methods, function(m) {
      draws <- theta_rec[[paste0("theta_", m)]]
      vals <- if (is_multivariate_theta) as.numeric(draws[, j]) else as.numeric(draws)
      data.frame(value = vals, method = method_labels[[m]])
    }))
    df$method <- factor(df$method, levels = unname(method_labels))
    p <- ggplot(df, aes(x = value, color = method)) +
      geom_density(linewidth = 0.8) +
      scale_color_manual(name = NULL, values = method_colors) +
      labs(x = theta_labels[[j]], y = NULL) +
      .sce3_theme()
    if (family == "bin") {
      p <- p + scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1))
    }
    if (!is.null(true_theta)) {
      true_theta_j <- if (is_multivariate_theta) true_theta[j] else true_theta[1]
      # grey, not black, so it doesn't visually collide with the (now
      # black) "NPP" curve -- same reasoning as elsewhere in this file.
      p <- p + geom_vline(xintercept = true_theta_j, linetype = "dashed",
                          color = "grey40", linewidth = 0.6)
    }
    p
  })
  theta_combined <- wrap_plots(theta_plots, ncol = n_theta)

  (eta_plot / theta_combined) +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom")
}

# =============================================================================
# General, family-agnostic boxplot of hat-eta / hat-theta across replicates
# =============================================================================
#
# Companion to plot_eta_theta_density(): that one shows the full posterior
# densities for a single replicate; this one shows -- averaged across every
# replicate in `sim` -- how hat(eta) and hat(theta) (the per-replicate
# posterior means, sim$hatdelta / sim$hattheta) are distributed for each
# method, i.e. "what happens on average". plot_sce_bin()/plot_sce_lm()
# already did this per-family; this version detects the family's structure
# automatically instead of hardcoding it, and reuses the existing
# plot_boxplot() helper -- which assumes the same column order sim_sce()
# always builds: hatdelta_npp/nppseq/onpp/onppseq, hattheta_npp/nppseq/
# onpp/onppseq (columns 1..4 = NPP, NPP-SEQ, ONPP, ONPP-SEQ).
#
# Arguments:
#   sim            one scenario's simulation object (same as
#                   plot_eta_theta_density(), e.g. sim_sces[[7]] for
#                   bin/normal, or sim_sces_beta_n0_ge_n[[2]] for lm)
#   family         one of "bin", "poi", "normal", "lm" -- only affects the
#                   title and how many true_theta reference lines to expect
#   true_theta     optional true theta value(s) to draw as a dotted
#                   reference line -- a single number for bin/poi/normal, or
#                   a vector matching hattheta's components for lm (e.g.
#                   c(betastar, sgstar))
#   scenario_label optional string (e.g. "Scenario III.II") shown in the title
#
# Returns a single patchwork object: one eta boxplot panel per component
# (hat(eta)[1], hat(eta)[2], ...) on top, theta boxplot panel(s) below --
# one per beta/sigma component for lm, a single panel otherwise -- each
# panel showing one box per method.
#
# Examples:
#   plot_eta_theta_boxplot(sim_sces[[7]], family = "bin", true_theta = 0.5,
#                           scenario_label = "Scenario III.II")
#   plot_eta_theta_boxplot(sim_sces_beta_n0_ge_n[[2]], family = "lm",
#                           true_theta = c(betastar, sgstar))
plot_eta_theta_boxplot <- function(sim, family = c("bin", "poi", "normal", "lm"),
                                    true_theta = NULL, scenario_label = NULL) {
  family <- match.arg(family)

  if (is.null(sim$hatdelta) || is.null(sim$hattheta)) {
    stop("`sim` doesn't look like a single scenario's simulation object -- ",
         "it has no $hatdelta/$hattheta. Did you pass the whole sim_sces ",
         "list instead of one element of it (e.g. sim_sces[[sce]])?")
  }

  # ---- eta boxplot panels: one per component -------------------------------
  n_eta <- length(sim$hatdelta)
  eta_plots <- lapply(seq_len(n_eta), function(i) {
    plot_boxplot(sim$hatdelta[[i]], bquote(hat(eta)[.(i)])) + coord_flip()
  })
  eta_combined <- wrap_plots(eta_plots, ncol = n_eta)

  # ---- theta boxplot panel(s) -----------------------------------------------
  # sim$hattheta is a single data.frame for bin/poi/normal (theta is a
  # scalar parameter) or a list of data.frames for lm (one per beta/sigma) --
  # detected here rather than hardcoded per family, same as
  # plot_eta_theta_density()'s is.matrix() check on the raw draws.
  is_multivariate_theta <- is.list(sim$hattheta) && !is.data.frame(sim$hattheta)
  n_theta <- if (is_multivariate_theta) length(sim$hattheta) else 1

  theta_labels <- if (!is_multivariate_theta) {
    list(expression(hat(theta)))
  } else {
    c(lapply(seq_len(n_theta - 1), function(i) bquote(hat(beta)[.(i)])),
      list(expression(hat(sigma))))
  }

  theta_plots <- lapply(seq_len(n_theta), function(j) {
    draws <- if (is_multivariate_theta) sim$hattheta[[j]] else sim$hattheta
    p <- plot_boxplot(draws, theta_labels[[j]]) + coord_flip()
    if (!is.null(true_theta)) {
      true_theta_j <- if (is_multivariate_theta) true_theta[j] else true_theta[1]
      p <- p + geom_hline(yintercept = true_theta_j, linetype = "dotted",
                          color = "red", linewidth = 0.6)
    }
    p
  })
  theta_combined <- wrap_plots(theta_plots, ncol = n_theta)

  title <- paste0(
    if (!is.null(scenario_label)) paste0(scenario_label, " — ") else "",
    "Average over replicates (", toupper(family), ")"
  )

  (eta_combined / theta_combined) +
    plot_annotation(title = title)
}

# =============================================================================
# Highlight only the best (smallest) value per row in a gt summary table
# =============================================================================
#
# Keeps the original per-row data_color() gradient (continuous mako/plasma
# coloring, anchored at each row's own min/max, same as before), with gt's
# automatic black-or-white per-cell text contrast left ON
# (autocolor_text = TRUE) -- the mako/plasma slice used here spans fairly
# dark to medium shades, and a single fixed text color was unreadable on
# the darker cells (see the "hard to read on dark background" report this
# function was patched for). On top of that automatic contrast, this
# function explicitly makes just the row's minimum bold + a distinct
# `best_text_color` (white by default), so exactly one number per
# row/group stands out as the best method, in addition to every other
# cell staying legible against its own background.
#
# Arguments:
#   tbl     a gt table object, already built (gt(), cols_label(), spanners,
#           etc. all applied)
#   df      the data.frame the table was built from -- row order must match
#           `tbl` (i.e. no filtering/reordering happened in between)
#   columns character vector of column names to color/highlight (typically
#           the method columns, e.g. c("NPP", "NPP_SEQ", "ONPP", "ONPP_SEQ"))
#   groups  optional list of character vectors partitioning `columns` into
#           groups that should each get their own row-wise gradient/best
#           cell, instead of one gradient across all of `columns`. Use this
#           when `columns` mixes quantities on different scales that
#           shouldn't share one color domain -- e.g. lm's summary tables
#           hold both beta and sigma columns per row, so beta and sigma
#           each need their own min/max: groups = list(beta_cols, sigma_cols).
#           Defaults to list(columns), i.e. one group covering everything.
#   palette     viridisLite palette function for the gradient (e.g.
#               viridisLite::mako for norm/MSE tables, viridisLite::plasma
#               for WIS tables -- matching the original)
#   best_text_color  text color forced onto the row's minimum cell
#
# Returns the gt table with the gradient and highlighting applied.
#
# Example:
#   table_mse_theta <- highlight_best_per_row(table_mse_theta, df_mse_theta,
#                                              columns = method_cols,
#                                              palette = viridisLite::mako)
#   table_wis_theta <- highlight_best_per_row(table_wis_theta, df_wis_theta,
#                                              columns = method_cols,
#                                              palette = viridisLite::plasma)
#   tbl <- highlight_best_per_row(tbl, df, columns = method_cols,
#                                  groups = list(beta_cols, sigma_cols),
#                                  palette = viridisLite::mako)
highlight_best_per_row <- function(tbl, df, columns, groups = NULL,
                                    palette = viridisLite::mako,
                                    best_text_color = "white",
                                    best_border_color = "black") {
  if (is.null(groups)) groups <- list(columns)

  # Relative luminance of a hex color (WCAG-style weights), used to pick
  # black vs white text ourselves instead of relying on gt's
  # autocolor_text: that heuristic turned out to flip inconsistently
  # between visually similar shades of the same gradient -- some cells
  # landed on black text on a background just as dark as a neighboring
  # cell that got white (the "nothing changed" report this was patched
  # for). Computing it ourselves from the exact rendered fill color below
  # is deterministic and doesn't depend on gt's internal contrast pick.
  .text_for_bg <- function(hex) {
    rgb <- grDevices::col2rgb(hex) / 255
    lin <- ifelse(rgb <= 0.04045, rgb / 12.92, ((rgb + 0.055) / 1.055)^2.4)
    L <- 0.2126 * lin[1, ] + 0.7152 * lin[2, ] + 0.0722 * lin[3, ]
    ifelse(L > 0.42, "black", "white")
  }

  for (i in seq_len(nrow(df))) {
    for (grp in groups) {
      row_vals <- as.numeric(df[i, grp])
      best_col <- grp[which.min(row_vals)]
      domain   <- c(min(row_vals), max(row_vals))
      pal      <- palette(100, direction = -1)[30:70]

      # Fill only -- no autocolor_text. Text color for every cell in this
      # group/row is set explicitly below, cell by cell.
      tbl <- tbl |>
        data_color(
          columns = all_of(grp),
          rows = i,
          palette = pal,
          domain = domain,
          autocolor_text = FALSE
        )

      # scales::col_numeric() with the same palette + domain reproduces
      # exactly the fill color data_color() just assigned to each cell
      # (gt uses this same scales function internally for a plain color
      # vector), so .text_for_bg() below is judging the actual rendered
      # background, not an approximation of it.
      color_fn <- scales::col_numeric(palette = pal, domain = domain)
      for (col in grp) {
        hex <- color_fn(as.numeric(df[i, col]))
        if (col == best_col) {
          # Best cell: forced best_text_color (bold) + border, so it
          # stands out regardless of how light/dark its own fill renders
          # -- same colorblind-safe discipline as the rest of this
          # project's palettes (see plot_theta_lm()'s comment in
          # analysis_seq_lm.r).
          style <- list(
            cell_text(color = best_text_color, weight = "bold"),
            cell_borders(sides = "all", color = best_border_color, weight = px(3))
          )
        } else {
          style <- cell_text(color = .text_for_bg(hex))
        }
        tbl <- tbl |>
          tab_style(style = style, locations = cells_body(columns = all_of(col), rows = i))
      }
    }
  }
  tbl
}

# =============================================================================
# Shared Monte Carlo replicates across scenarios that share a theta0 *set*
# =============================================================================
#
# Motivating case: scenarios_bin.R's "cancellation effect" pair,
#   "5_1" = c(0.10, 0.90, 0.50), "5_2" = c(0.50, 0.90, 0.10)
# is the same three theta0 values with positions 1 and 3 swapped. NPP and
# NPP-SEQ's eta prior (eta_bin.stan) is iid/exchangeable across the K
# historical datasets, with no ordering constraint -- so relabeling which
# position holds which theta0 value is a pure symmetry for them: their
# *true* MSE/WIS must be identical between "5_1" and "5_2". ONPP/ONPP-SEQ
# (gamma_bin.stan's eta = cumulative_sum(gamma), forcing eta1<=eta2<=eta3)
# are genuinely order-sensitive and are NOT expected to match -- that's the
# whole point of testing this permutation.
#
# Today each scenario draws its own historical data (z0) and current data
# (z) independently inside sample_sce_bin(), even when seeded identically,
# because data is drawn *by position* (position k's rbinom() call always
# happens before position k+1's, consuming a different point in the RNG
# stream), not by theta0 *value*. So even with a shared seed, "5_1"'s
# position-1 draw (theta0 = 0.10) and "5_2"'s position-3 draw (also
# theta0 = 0.10) land at different stream offsets and get different
# numbers -- same seed, but not the same data. That leaves NPP/NPP-SEQ's
# reported MSE/WIS differing by ordinary Monte Carlo noise between "5_1"
# and "5_2", instead of being (near-)identical the way the underlying
# symmetry says they should be.
#
# build_shared_replicate_data() fixes this at the data level: for each
# *unique* theta0 value appearing anywhere in the group, it draws num_sim
# iid replicates ONCE, then reassembles each scenario's per-replicate z0 by
# looking up the shared draw for the value at each of ITS OWN positions --
# so scenario "5_1"'s position 1 (theta0 = 0.10) and scenario "5_2"'s
# position 3 (theta0 = 0.10) get literally the same draw for a given
# replicate index, correctly relabeled rather than reused positionally.
# The current-study data z (depends only on n/theta, not theta0) is drawn
# once and shared identically across every scenario in the group.
#
# Feed the result into sim_sce()'s `replicate_data` argument (one call per
# scenario, using that scenario's own entry) to make its data generation
# use these shared draws instead of drawing fresh internally.
#
# NOTE on what this does and doesn't guarantee: with shared data, NPP/
# NPP-SEQ's *posterior given the data* is exactly the same distribution
# for matched replicates across the group (the data is identical up to
# relabeling, and the model is exchangeable) -- but .fit_family() doesn't
# pass an explicit seed to cmdstanr's $sample() calls, so each scenario's
# Stan run still draws its own independent MCMC seed. That means results
# won't be bit-for-bit identical even after this fix, but the dominant
# source of divergence (independently-simulated data) is eliminated,
# leaving only ordinary MCMC sampling noise -- which should be much
# smaller. Forcing fully deterministic, bit-identical output would
# additionally require passing a shared seed into $sample() itself, which
# this function does not do (that's a separate, bigger change to
# .fit_family()/sample_sce_bin(), not made here).
#
# Arguments:
#   theta0_list  named list, one theta0 vector per scenario in the group
#                (e.g. list(`5_1` = c(0.10, 0.90, 0.50),
#                           `5_2` = c(0.50, 0.90, 0.10))) -- every vector
#                must be a permutation of the same multiset of values, or
#                this stops with an error (sharing draws only makes sense
#                across scenarios that actually share a theta0 set).
#   n0           the (shared) per-position sample sizes, e.g. c(100,100,100)
#   n, theta     the current-study sample size and true theta -- must be the
#                same for every scenario in the group (asserted implicitly:
#                pass the single shared values, not a per-scenario list)
#   num_sim      number of replicates to generate
#   seed         optional seed for reproducibility (set.seed() on the
#                current process -- this runs sequentially, not on a
#                cluster, so a plain set.seed() suffices here)
#
# Returns a named list (same names as theta0_list), each entry a list of
# length num_sim of list(z0 = <vector, length(n0)>, z = <scalar>), ready to
# pass as sim_sce()'s `replicate_data` for that scenario.
#
# Example:
#   shared_v <- build_shared_replicate_data(
#     theta0_list = list(`5_1` = scenarios[["5_1"]]$theta0,
#                        `5_2` = scenarios[["5_2"]]$theta0),
#     n0 = scenarios[["5_1"]]$n0, n = scenarios[["5_1"]]$n,
#     theta = scenarios[["5_1"]]$theta,
#     num_sim = n_sim, seed = base_seed
#   )
#   sims[["5_1"]] <- sim_sce(model = "bin", n_cores, n_sim, scenarios[["5_1"]],
#                             gamma_model_bin, delta_model_bin,
#                             replicate_data = shared_v[["5_1"]])
#   sims[["5_2"]] <- sim_sce(model = "bin", n_cores, n_sim, scenarios[["5_2"]],
#                             gamma_model_bin, delta_model_bin,
#                             replicate_data = shared_v[["5_2"]])
build_shared_replicate_data <- function(theta0_list, n0, n, theta, num_sim, seed = NULL) {
  scenario_names <- names(theta0_list)
  if (is.null(scenario_names) || any(scenario_names == "")) {
    stop("theta0_list must be a named list, one theta0 vector per scenario.")
  }
  K <- length(n0)
  for (nm in scenario_names) {
    if (length(theta0_list[[nm]]) != K) {
      stop(sprintf(
        "theta0_list[['%s']] has length %d but n0 has length %d -- every scenario's ",
        "theta0 must have one entry per historical dataset in n0.",
        nm, length(theta0_list[[nm]])
      ))
    }
  }

  # Every scenario's theta0 must be a permutation of the SAME multiset of
  # values -- otherwise "shared replicates" doesn't mean anything, and
  # silently proceeding would quietly share draws between scenarios that
  # don't actually agree on their theta0 set.
  ref_sorted <- sort(theta0_list[[1]])
  for (nm in scenario_names[-1]) {
    if (!isTRUE(all.equal(sort(theta0_list[[nm]]), ref_sorted))) {
      stop(sprintf(
        "theta0_list[['%s']] = (%s) is not a permutation of theta0_list[['%s']] = (%s) -- ",
        "build_shared_replicate_data() only makes sense for scenarios that share the ",
        "exact same set of theta0 values, just assigned to different positions.",
        nm, paste(theta0_list[[nm]], collapse = ", "),
        scenario_names[1], paste(ref_sorted, collapse = ", ")
      ))
    }
  }

  # Every occurrence of a given theta0 value must be paired with the same
  # n0, or the "shared draw" for that value wouldn't actually be drawn
  # from the same distribution at every position it appears in. With the
  # current n0 = c(100,100,100) shared across every position this always
  # holds, but this guards against a future scenario table where n0 varies
  # by position.
  n0_for_value <- list()
  for (nm in scenario_names) {
    theta0 <- theta0_list[[nm]]
    for (k in seq_len(K)) {
      key <- as.character(theta0[k])
      if (is.null(n0_for_value[[key]])) {
        n0_for_value[[key]] <- n0[k]
      } else if (n0_for_value[[key]] != n0[k]) {
        stop(sprintf(
          "theta0 value %s has n0 = %d at one position but n0 = %d at another ",
          "(scenario '%s', position %d) -- shared replicates require the same ",
          "n0 wherever a given theta0 value appears.",
          key, n0_for_value[[key]], n0[k], nm, k
        ))
      }
    }
  }

  if (!is.null(seed)) set.seed(seed)

  # One shared draw per unique theta0 value, generated once and reused --
  # relabeled into position -- by every scenario that needs it.
  z0_by_value <- lapply(names(n0_for_value), function(key) {
    rbinom(num_sim, n0_for_value[[key]], as.numeric(key))
  })
  names(z0_by_value) <- names(n0_for_value)

  # Current-study data doesn't depend on theta0 at all -- same n/theta for
  # every scenario in the group, so it's drawn once and shared identically.
  z_shared <- rbinom(num_sim, n, theta)

  # Assemble each scenario's per-replicate z0 by VALUE, not position: this
  # is the step that has to get the order right. Scenario nm's position k
  # gets the shared draw for theta0_list[[nm]][k] -- whatever value that
  # scenario actually has there -- never "whatever happened to be at
  # position k" in some other scenario.
  result <- lapply(scenario_names, function(nm) {
    theta0 <- theta0_list[[nm]]
    lapply(seq_len(num_sim), function(i) {
      z0_i <- vapply(seq_len(K), function(k) z0_by_value[[as.character(theta0[k])]][i], numeric(1))
      list(z0 = z0_i, z = z_shared[i])
    })
  })
  names(result) <- scenario_names
  result
}
