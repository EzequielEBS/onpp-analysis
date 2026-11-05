get_summary_lm <- function(sim_sces, true_value) {
  # NOTE: this used to also load("../../simulated-data-regression/data/
  # true_params.RData") here, but nothing in this function ever referenced
  # any object from that file -- it only uses the `true_value` argument the
  # caller passes in. Removed; callers that need the true parameters (to
  # build `true_value`) already load that file themselves (see
  # analysis_lm.R/analysis_seq_lm.r), so this was a pure no-op re-read of
  # the same file on every call.
  #
  # `true_value` is c(beta1, beta2, beta3, sigma): the first 3 entries of
  # `sim$hattheta` are always the regression coefficients and the last is
  # sigma (see analysis_lm.R, `true_value = c(betastar, sgstar)`).
  #
  # MSE instead of norm: for beta (3 coefficients) this is the per-component
  # mean squared error, averaged over the 3 coefficients (not the
  # Euclidean norm of the error vector); for sigma (a scalar) it's just
  # the mean squared error of that one component, the 1-dimensional
  # analogue.
  mses <- lapply(sim_sces, function(sim) {
    n_params <- length(sim$hattheta)
    beta_idx  <- seq_len(n_params - 1)
    sigma_idx <- n_params

    errors <- lapply(seq_len(n_params), function(i) {
      as.matrix(sim$hattheta[[i]]) - true_value[i]
    })

    mse_beta  <- Reduce(`+`, lapply(errors[beta_idx], function(e) colMeans(e^2))) / length(beta_idx)
    mse_sigma <- colMeans(errors[[sigma_idx]]^2)

    list(mse_beta = mse_beta, mse_sigma = mse_sigma)
  })

  biass <- lapply(sim_sces, function(sim) {
    do.call(rbind, lapply(seq_along(sim$hattheta), function(i) {
      colMeans((sim$hattheta[[i]] - true_value[i]))
    }))
  })

  quantile_level <- c(0.025, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.975)
  wis_list <- lapply(sim_sces, function(sim) {
    compute_wis_lm(sim, quantile_level, true_value)
  })

  # ------------------------------------------------------------------------------
  # Summary results
  # ------------------------------------------------------------------------------
  #
  # WIS_beta below is the mean of the per-coordinate WIS over beta1/beta2/
  # beta3 (indices 1:3 of each hat_wis_* vector); WIS_sigma is the WIS for
  # sigma alone (index 4) -- mirroring the beta/sigma split used for the
  # norm above, instead of the single all-4-parameters-combined WIS this
  # used to report.

  summary_table <- data.frame(scenario = unlist(lapply(1:length(sim_sces),
                                                      function(i) {
                                                        sce <- as.roman(ceiling(i/3))
                                                        sce <- ifelse(i %% 3 == 1, paste0(sce,".I"), ifelse(i %% 3 == 2, paste0(sce,".II"), paste0(sce,".III")))
                                                        rep(sce,4)
                                                      })),
                              model = rep(c("NPP", "NPP-SEQ", "ONPP", "ONPP-SEQ"),length(sim_sces)),
                              MSE_beta = unlist(lapply(mses, function(x) x$mse_beta)),
                              MSE_sigma = unlist(lapply(mses, function(x) x$mse_sigma)),
                              bias_beta1 = unlist(lapply(biass, function(x) x[1,])),
                              bias_beta2 = unlist(lapply(biass, function(x) x[2,])),
                              bias_beta3 = unlist(lapply(biass, function(x) x[3,])),
                              bias_sg = unlist(lapply(biass, function(x) x[4,])),
                              WIS_beta = unlist(lapply(wis_list, function(x) c(mean(x$hat_wis_npp[1:3]),
                                                                                mean(x$hat_wis_nppseq[1:3]),
                                                                                mean(x$hat_wis_onpp[1:3]),
                                                                                mean(x$hat_wis_onppseq[1:3])))),
                              WIS_sigma = unlist(lapply(wis_list, function(x) c(x$hat_wis_npp[4],
                                                                                x$hat_wis_nppseq[4],
                                                                                x$hat_wis_onpp[4],
                                                                                x$hat_wis_onppseq[4])))
  )

  row.names(summary_table) <- NULL

  print(xtable::xtable(summary_table %>%
                        mutate(
                          MSE_beta = formatC(MSE_beta, format = "f", digits = 3),
                          MSE_sigma = formatC(MSE_sigma, format = "f", digits = 3),
                          bias_beta1 = formatC(bias_beta1, format = "f", digits = 3),
                          bias_beta2 = formatC(bias_beta2, format = "f", digits = 3),
                          bias_beta3 = formatC(bias_beta3, format = "f", digits = 3),
                          bias_sg = formatC(bias_sg, format = "f", digits = 3),
                          WIS_beta = formatC(WIS_beta, format = "f", digits = 3),
                          WIS_sigma = formatC(WIS_sigma, format = "f", digits = 3)
                        )
  ),
  type = "latex", include.rownames = FALSE)
  return(list(summary_table = summary_table, mses = mses, biass = biass, wis_list = wis_list))
}
