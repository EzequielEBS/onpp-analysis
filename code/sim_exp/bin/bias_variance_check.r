# Bias/variance/MSE decomposition for scenario 13 (V.I), verified against
# get_summary_bin.r's own mses/biass computation.
#
# Run this after sim_sces / true_value are loaded (i.e. after sourcing
# analysis_bin.r up through the point where sim_sces is built).

sim13 <- sim_sces[[14]]
methods <- c("NPP", "NPP-SEQ", "ONPP", "ONPP-SEQ")

hattheta <- sim13$hattheta  # data.frame: hattheta_npp/nppseq/onpp/onppseq

bias <- colMeans(hattheta) - true_value
variance <- apply(hattheta, 2, var)
# population-style MSE (matches get_summary_bin.r's colMeans((hattheta -
# true_value)^2), which divides by n, not n-1 -- so we reconstruct it the
# same way here for an apples-to-apples comparison, and separately via the
# bias^2 + variance identity using the n-1 sample variance from var(), to
# show both land in the same place)
mse_direct <- colMeans((as.matrix(hattheta) - true_value)^2)
n <- nrow(hattheta)
mse_from_decomp <- bias^2 + variance * (n - 1) / n  # matches population variance

decomp <- data.frame(
  method = methods,
  bias = round(as.numeric(bias), 4),
  variance = round(as.numeric(variance), 5),
  mse_direct = round(as.numeric(mse_direct), 5),
  mse_from_bias2_plus_var = round(as.numeric(mse_from_decomp), 5)
)
print(decomp)


cat("\nOutlier check -- how far the single worst replicate is from the truth,\n",
    "per method (this is what the boxplot's dots show and MSE weights\n",
    "quadratically, but a box's IQR doesn't capture at all):\n", sep = "")
worst <- sapply(hattheta, function(x) max(abs(x - true_value)))
print(round(worst, 4))

cat("\nCross-check against get_summary_bin()'s own numbers for this scenario\n",
    "(requires having already run summary_bin <- get_summary_bin(sim_sces)):\n", sep = "")
if (exists("summary_bin")) {
  print(summary_bin$summary_table[summary_bin$summary_table$Scenario == "V.I", ])
} else {
  cat("(summary_bin not found in this session -- run summary_bin <- get_summary_bin(sim_sces) first if you want this cross-check)\n")
}
