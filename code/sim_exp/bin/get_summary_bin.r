get_summary_bin <- function(sim_sces, true_value = 0.5) {
  # Compute MSE
  mses <- lapply(sim_sces, function(sim) {
    colMeans((sim$hattheta - true_value)^2)
  })

  biass <- lapply(sim_sces, function(sim) {
    colMeans(sim$hattheta) - true_value
  })


  # Compute BCI
  bcis_95 <- lapply(sim_sces, function(sim) {
    compute_bci(sim$theta, 0.95)
  })
  bcis_90 <- lapply(sim_sces, function(sim) {
    compute_bci(sim$theta, 0.90)
  })
  bcis_80 <- lapply(sim_sces, function(sim) {
    compute_bci(sim$theta, 0.80)
  })
  bcis_50 <- lapply(sim_sces, function(sim) {
    compute_bci(sim$theta, 0.50)
  })
  # compute coverage
  coverages_95 <- lapply(1:length(bcis_95), function(i) {
    compute_coverage(bcis_95[[i]], true_value)
  })
  coverages_90 <- lapply(1:length(bcis_90), function(i) {
    compute_coverage(bcis_90[[i]], true_value)
  })
  coverages_80 <- lapply(1:length(bcis_80), function(i) {
    compute_coverage(bcis_80[[i]], true_value)
  })
  coverages_50 <- lapply(1:length(bcis_50), function(i) {
    compute_coverage(bcis_50[[i]], true_value)
  })

  # compute WIS
  quantile_level <- c(0.025, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.975)
  wis_list <- lapply(sim_sces, function(sim) {
    compute_wis(sim$theta, quantile_level, true_value)
  })

  # ------------------------------------------------------------------------------
  # Summary results
  # ------------------------------------------------------------------------------

  summary_table <- data.frame(scenario = unlist(lapply(1:length(sim_sces), 
                                                      function(i) {
                    sce <- as.roman(ceiling(i/3))
                    sce <- ifelse(i %% 3 == 1, paste0(sce,".I"), ifelse(i %% 3 == 2, paste0(sce,".II"), paste0(sce,".III")))
                    rep(sce,4)
                  })
                  ),
                  model = rep(c("NPP", "NPP-SEQ", "ONPP", "ONPP-SEQ"),length(sim_sces)),
                  bias = unlist(biass),
                  mse = unlist(mses),
                  average_wis = unlist(wis_list),
                  average_len_95 = unlist(lapply(bcis_95, compute_avg_len)),
                  average_len_90 = unlist(lapply(bcis_90, compute_avg_len)),
                  average_len_80 = unlist(lapply(bcis_80, compute_avg_len)),
                  average_len_50 = unlist(lapply(bcis_50, compute_avg_len)),
                  coverage_95 = unlist(coverages_95),
                  coverage_90 = unlist(coverages_90),
                  coverage_80 = unlist(coverages_80),
                  coverage_50 = unlist(coverages_50)
  )
  colnames(summary_table) <- c("Scenario", "Prior", "Bias", "MSE", "Avg_WIS",
                                "Avg_Len_BCI95", "Avg_Len_BCI90", "Avg_Len_BCI80", "Avg_Len_BCI50",
                              "Cov_BCI95", "Cov_BCI90", "Cov_BCI80", "Cov_BCI50")
  # NOTE: this used to unconditionally blank the roman-numeral scenario
  # labels computed just above (`summary_table$Scenario <- ""`), which threw
  # away real information for no apparent reason. Removed so the printed
  # xtable actually shows which scenario each row belongs to.
  row.names(summary_table) <- NULL

  print(xtable::xtable(summary_table %>% 
                        mutate(
                          MSE = formatC(MSE, format = "f", digits = 3),
                          Bias = formatC(Bias, format = "f", digits = 3),
                          Avg_WIS = formatC(Avg_WIS, format = "f", digits = 3)
                        )
  ), 
  type = "latex", include.rownames = FALSE)
  return(list(
    summary_table = summary_table, mses = mses, biass = biass, wis_list = wis_list,
    bcis_95 = bcis_95, bcis_90 = bcis_90, bcis_80 = bcis_80, bcis_50 = bcis_50,
    coverages_95 = coverages_95, coverages_90 = coverages_90, coverages_80 = coverages_80, 
    coverages_50 = coverages_50
  ))
}