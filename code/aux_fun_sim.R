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

  theta_draws <- list(
    npp     = rbeta(nrow(eta_draws$npp),
                     z + a + eta_draws$npp %*% z0,
                     n - z + b + eta_draws$npp %*% (n0 - z0)),
    nppseq  = rbeta(nrow(eta_draws$nppseq),
                     z + K * (a - 1) + 1 + eta_draws$nppseq %*% z0,
                     n - z + K * (b - 1) + 1 + eta_draws$nppseq %*% (n0 - z0)),
    onpp    = rbeta(nrow(eta_draws$onpp),
                     z + a + eta_draws$onpp %*% z0,
                     n - z + b + eta_draws$onpp %*% (n0 - z0)),
    onppseq = rbeta(nrow(eta_draws$onppseq),
                     z + K * (a - 1) + 1 + eta_draws$onppseq %*% z0,
                     n - z + K * (b - 1) + 1 + eta_draws$onppseq %*% (n0 - z0))
  )

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
sim_sce <- function(model, num_cores, num_sim, sce, gamma_model, delta_model, seed = NULL) {
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
    results <- pbapply::pblapply(
      seq_len(num_sim),
      function(i, sce, gamma_model, delta_model) {
        tryCatch(
          sample_sce(sce, gamma_model = gamma_model, delta_model = delta_model),
          error = function(e) list(.failed = TRUE, .error = conditionMessage(e))
        )
      },
      sce = sce,
      gamma_model = gamma_model,
      delta_model = delta_model,
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

plot_sce3_theta <- function(sub_sce, sim, sim_sces, true_value = get("true_value", envir = .GlobalEnv)) {
  # true_value defaults to looking up the global of the same name (see the
  # note in plot_sce_bin() above) so existing call sites keep working.
  name_sce <- as.roman(ceiling(3))
  name_sce <- ifelse(sub_sce %% 3 == 1, paste0(name_sce,".I"),
                     ifelse(sub_sce %% 3 == 2, paste0(name_sce,".II"), paste0(name_sce,".III")))
  sim_sub_sce <- sim_sces[[sub_sce]]
  plot_theta <- lapply(1:length(sim_sub_sce), function(i) {
    plot <- ggplot() +
    geom_density(aes(x = sim_sub_sce[[i]]$theta[[sim]]$theta_npp, color = "NPP"), linewidth = 0.8) +
    geom_density(aes(x = sim_sub_sce[[i]]$theta[[sim]]$theta_onpp, color = "ONPP"), linewidth = 0.8) +
    geom_vline(xintercept = true_value, linetype = "dotted", color = "black", size = 1) +
    scale_color_manual(name = NULL,
                       values = c("NPP" = RColorBrewer::brewer.pal(4, "Set2")[1],
                                  "ONPP" = RColorBrewer::brewer.pal(4, "Set2")[2])) +
    labs(x = expression(theta),
         y = NULL) +
    xlim(0,1)
    if (i == 1) {
      plot <- plot +
          ggtitle(label = paste0("Scenario ", name_sce),
                  subtitle = paste0("(", LETTERS[i], ")"))
    } else {
      plot <- plot +
        ggtitle(label = NULL,
                subtitle = paste0("(", LETTERS[i], ")"))
    }
    return(plot)
  })
  combined_plots <- wrap_plots(plot_theta, ncol = length(plot_theta)) +
    plot_layout(guides = "collect")
  return(combined_plots)
  }

plot_sce3_eta <- function(sub_sce, sim, sim_sces){
  name_sce <- as.roman(ceiling(3))
  name_sce <- ifelse(sub_sce %% 3 == 1, paste0(name_sce,".I"), 
                     ifelse(sub_sce %% 3 == 2, paste0(name_sce,".II"), paste0(name_sce,".III")))
  
  sim_sub_sce <- sim_sces[[sub_sce]]
  plots_npp <- lapply(1:length(sim_sub_sce), function(i) {
    plot <- ggplot() +
      geom_density(bounds = c(0, 1), aes(x = sim_sub_sce[[i]]$delta[[sim]]$delta_npp[,1], color = "delta1"), linewidth = 0.8) +
      geom_density(bounds = c(0, 1), aes(x = sim_sub_sce[[i]]$delta[[sim]]$delta_npp[,2], color = "delta2"), linewidth = 0.8) +
      geom_density(bounds = c(0, 1), aes(x = sim_sub_sce[[i]]$delta[[sim]]$delta_npp[,3], color = "delta3"), linewidth = 0.8) +
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
          ggtitle(label = paste0("Scenario ", name_sce, "\nNPP"),
                  subtitle = paste0("(", LETTERS[i], ")"))
    } else {
      plot <- plot + 
        ggtitle(label = NULL,
                subtitle = paste0("(", LETTERS[i], ")"))
    }
    return(plot)
    })
  
  plots_onpp <- lapply(1:length(sim_sub_sce), function(i) {
    plot <- ggplot() +
      geom_density(bounds = c(0, 1), aes(x = sim_sub_sce[[i]]$delta[[sim]]$delta_onpp[,1], color = "delta1"), linewidth = 0.8) +
      geom_density(bounds = c(0, 1), aes(x = sim_sub_sce[[i]]$delta[[sim]]$delta_onpp[,2], color = "delta2"), linewidth = 0.8) +
      geom_density(bounds = c(0, 1), aes(x = sim_sub_sce[[i]]$delta[[sim]]$delta_onpp[,3], color = "delta3"), linewidth = 0.8) +
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
  
  combined_plots <- wrap_plots(plots_npp, ncol = length(plots_npp)) / 
    wrap_plots(plots_onpp, ncol = length(plots_onpp)) +
    plot_layout(guides = "collect")
  return(combined_plots)
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
