# Scenario definitions for the binomial NPP/ONPP simulation study
# (code/sim_exp/bin/). This is the single source of truth for every
# scenario's prior hyperparameters (n0, n, a, b, al, bl, alpha, theta0) and
# for the output file paths those scenarios are saved to/read from.
#
# It is sourced by two scripts that used to define these tables
# independently (and drift out of sync):
#   - simulations_bin.R runs sim_sce()/sample_sce_bin() over every entry in
#     `scenarios` and saves results to `scenario_files[[key]]`.
#   - analysis_bin.r reads those same results back and, for a couple of
#     closed-form prior plots, needs the scenario's own $a/$b hyperparameters
#     (e.g. scenarios[["3_1_all_flat"]]$a) rather than re-deriving them by
#     hand as untracked globals.
#
# Naming follows the paper's scenario numbering: I = compatible order,
# II = almost compatible, III = incompatible, IV = neutral,
# V = cancellation effect; .1/.2/.3 = high/small/no congruence (V only has
# two sub-scenarios).

base_par <- list(n0 = c(100, 100, 100), n = 100,
                  a = 1/2, b = 1/2, al = 1/2, bl = 1/2,
                  alpha = rep(1/4, 4))

# Main scenarios: default priors, all 14 (order x congruence) settings.
main_theta0 <- list(
  "1_1" = c(0.30, 0.40, 0.50), "1_2" = c(0.20, 0.30, 0.40), "1_3" = c(0.05, 0.15, 0.25),
  "2_1" = c(0.30, 0.50, 0.40), "2_2" = c(0.20, 0.40, 0.30), "2_3" = c(0.05, 0.25, 0.15),
  "3_1" = c(0.50, 0.40, 0.30), "3_2" = c(0.40, 0.30, 0.20), "3_3" = c(0.25, 0.15, 0.05),
  "4_1" = c(0.50, 0.50, 0.50), "4_2" = c(0.30, 0.30, 0.30), "4_3" = c(0.10, 0.10, 0.10),
  "5_1" = c(0.10, 0.90, 0.50), "5_2" = c(0.50, 0.90, 0.10)
)

# Scenario III (incompatible order) is re-run under three progressively
# flatter prior specifications, to see how much the eta/gamma priors matter,
# using the same three theta0 settings as sce3.1/3.2/3.3 above.
sce3_theta0 <- main_theta0[c("3_1", "3_2", "3_3")]
flat_variants <- list(
  eta_flat       = list(overrides = list(al = 1, bl = 1),
                        suffix = "_eta_flat"),
  eta_gamma_flat = list(overrides = list(al = 1, bl = 1, alpha = rep(1, 4)),
                        suffix = "_eta_gamma_flat"),
  all_flat       = list(overrides = list(a = 1, b = 1, al = 1, bl = 1, alpha = rep(1, 4)),
                        suffix = "_all_flat")
)

scenarios <- list()
scenario_files <- character(0)

for (id in names(main_theta0)) {
  scenarios[[id]] <- modifyList(base_par, list(theta0 = main_theta0[[id]], theta = 0.5))
  scenario_files[id] <- paste0("results/samples/bin/sim", id, "_bin.qs2")
}

for (v in names(flat_variants)) {
  variant <- flat_variants[[v]]
  for (id in names(sce3_theta0)) {
    key <- paste0(id, "_", v)
    scenarios[[key]] <- modifyList(modifyList(base_par, variant$overrides),
                                    list(theta0 = sce3_theta0[[id]], theta = 0.5))
    scenario_files[key] <- paste0("results/samples/bin/sim", id, "_bin", variant$suffix, ".qs2")
  }
}
