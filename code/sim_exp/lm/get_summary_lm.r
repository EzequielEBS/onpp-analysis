get_summary_lm <- function(sim_sces, true_value) {
  # NOTE: this used to also load("../../simulated-data-regression/data/
  # true_params.RData") here, but nothing in this function ever referenced
  # any object from that file -- it only uses the `true_value` argument the
  # caller passes in. Removed; callers that need the true parameters (to
  # build `true_value`) already load that file themselves (see
  # analysis_lm.R/analysis_seq_lm.r), so this was a pure no-op re-read of
  # the same file on every call.
  #
  # This used to compute the same per-model MSE vector twice, in two
  # different orders (mean-over-replicates-of-sum-over-params vs.
  # sum-over-params-of-mean-over-replicates-per-param -- mathematically
  # identical by linearity, verified numerically before removing). The
  # first computation's result was never assigned to anything and was
  # immediately discarded; only the second (kept below) was ever returned.
  mses <- lapply(sim_sces, function(sim) {
    colSums(
      do.call(rbind, lapply(seq_along(sim$hattheta), function(i) {
        colMeans((sim$hattheta[[i]] - true_value[i])^2)
      }))
    )
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

  summary_table <- data.frame(scenario = unlist(lapply(1:length(sim_sces), 
                                                      function(i) {
                                                        sce <- as.roman(ceiling(i/3))
                                                        sce <- ifelse(i %% 3 == 1, paste0(sce,".I"), ifelse(i %% 3 == 2, paste0(sce,".II"), paste0(sce,".III")))
                                                        rep(sce,4)
                                                      })),
                              model = rep(c("NPP", "NPP-SEQ", "ONPP", "ONPP-SEQ"),length(sim_sces)),
                              MSE = unlist(mses),
                              bias_beta1 = unlist(lapply(biass, function(x) x[1,])),
                              bias_beta2 = unlist(lapply(biass, function(x) x[2,])),
                              bias_beta3 = unlist(lapply(biass, function(x) x[3,])),
                              bias_sg = unlist(lapply(biass, function(x) x[4,])),
                              WIS_beta1 = unlist(lapply(wis_list, function(x) c(x$hat_wis_npp[1],
                                                                                x$hat_wis_nppseq[1],
                                                                                x$hat_wis_onpp[1],
                                                                                x$hat_wis_onppseq[1]))),
                              WIS_beta2 = unlist(lapply(wis_list, function(x) c(x$hat_wis_npp[2],
                                                                                x$hat_wis_nppseq[2],
                                                                                x$hat_wis_onpp[2],
                                                                                x$hat_wis_onppseq[2]))),
                              WIS_beta3 = unlist(lapply(wis_list, function(x) c(x$hat_wis_npp[3],
                                                                                x$hat_wis_nppseq[3],
                                                                                x$hat_wis_onpp[3],
                                                                                x$hat_wis_onppseq[3]))),
                              WIS_sg = unlist(lapply(wis_list, function(x) c(x$hat_wis_npp[4],
                                                                              x$hat_wis_nppseq[4],
                                                                              x$hat_wis_onpp[4],
                                                                              x$hat_wis_onppseq[4])))
  )

  row.names(summary_table) <- NULL

  print(xtable::xtable(summary_table %>% 
                        mutate(
                          MSE = formatC(MSE, format = "f", digits = 3),
                          bias_beta1 = formatC(bias_beta1, format = "f", digits = 3),
                          bias_beta2 = formatC(bias_beta2, format = "f", digits = 3),
                          bias_beta3 = formatC(bias_beta3, format = "f", digits = 3),
                          bias_sg = formatC(bias_sg, format = "f", digits = 3),
                          WIS_beta1 = formatC(WIS_beta1, format = "f", digits = 3),
                          WIS_beta2 = formatC(WIS_beta2, format = "f", digits = 3),
                          WIS_beta3 = formatC(WIS_beta3, format = "f", digits = 3),
                          WIS_sg = formatC(WIS_sg, format = "f", digits = 3)
                        )
  ), 
  type = "latex", include.rownames = FALSE)
  return(list(summary_table = summary_table, mses = mses, biass = biass, wis_list = wis_list))
}
