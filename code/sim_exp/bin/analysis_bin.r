# Load libraries
library(cmdstanr)
library(bayesplot)
library(MCMCpack)
library(dplyr)
library(tidyverse)
library(posterior)
library(patchwork)
library(parallel)
library(latex2exp)
library(scoringutils)
library(gt)
library(gtExtras)
library(ggh4x)  # needed for plot_sce3_eta()'s facet_grid2(independent = "y")

# run auxiliary functions
source("code/aux_fun_sim.R")
source("code/sim_exp/bin/get_summary_bin.r")
# scenario definitions (n0/n/a/b/al/bl/alpha/theta0 per scenario) -- shared
# with simulations_bin.R; used below to recover e.g. scenarios[["3_1_all_flat"]]$a
# for the closed-form theta prior plot, instead of relying on untracked globals.
source("code/sim_exp/bin/scenarios_bin.R")

# load sample
sim1.1 <- qs2::qs_read("results/samples/bin/sim1_1_bin.qs2")
sim1.2 <- qs2::qs_read("results/samples/bin/sim1_2_bin.qs2")
sim1.3 <- qs2::qs_read("results/samples/bin/sim1_3_bin.qs2")
sim2.1 <- qs2::qs_read("results/samples/bin/sim2_1_bin.qs2")
sim2.2 <- qs2::qs_read("results/samples/bin/sim2_2_bin.qs2")
sim2.3 <- qs2::qs_read("results/samples/bin/sim2_3_bin.qs2")
sim3.1 <- qs2::qs_read("results/samples/bin/sim3_1_bin.qs2")
sim3.2 <- qs2::qs_read("results/samples/bin/sim3_2_bin.qs2")
sim3.3 <- qs2::qs_read("results/samples/bin/sim3_3_bin.qs2")
sim4.1 <- qs2::qs_read("results/samples/bin/sim4_1_bin.qs2")
sim4.2 <- qs2::qs_read("results/samples/bin/sim4_2_bin.qs2")
sim4.3 <- qs2::qs_read("results/samples/bin/sim4_3_bin.qs2")
sim5.1 <- qs2::qs_read("results/samples/bin/sim5_1_bin.qs2")
sim5.2 <- qs2::qs_read("results/samples/bin/sim5_2_bin.qs2")

sim3.1_eta_flat <- qs2::qs_read("results/samples/bin/sim3_1_bin_eta_flat.qs2")
sim3.2_eta_flat <- qs2::qs_read("results/samples/bin/sim3_2_bin_eta_flat.qs2")
sim3.3_eta_flat <- qs2::qs_read("results/samples/bin/sim3_3_bin_eta_flat.qs2")

sim3.1_eta_gamma_flat <- qs2::qs_read("results/samples/bin/sim3_1_bin_eta_gamma_flat.qs2")
sim3.2_eta_gamma_flat <- qs2::qs_read("results/samples/bin/sim3_2_bin_eta_gamma_flat.qs2")
sim3.3_eta_gamma_flat <- qs2::qs_read("results/samples/bin/sim3_3_bin_eta_gamma_flat.qs2")

sim3.1_all_flat <- qs2::qs_read("results/samples/bin/sim3_1_bin_all_flat.qs2")
sim3.2_all_flat <- qs2::qs_read("results/samples/bin/sim3_2_bin_all_flat.qs2")
sim3.3_all_flat <- qs2::qs_read("results/samples/bin/sim3_3_bin_all_flat.qs2")

samples_prior_sce3.1_all_flat <- qs2::qs_read("results/samples/bin/samples_prior_sce3_1_all_flat.qs2")
samples_prior_sce3.2_all_flat <- qs2::qs_read("results/samples/bin/samples_prior_sce3_2_all_flat.qs2")
samples_prior_sce3.3_all_flat <- qs2::qs_read("results/samples/bin/samples_prior_sce3_3_all_flat.qs2")

# true value
true_value <- 0.5

# combine simulations
sim_sces <- list(
  sim1.1, sim1.2, sim1.3,
  sim2.1, sim2.2, sim2.3,
  sim3.1, sim3.2, sim3.3,
  sim4.1, sim4.2, sim4.3,
  sim5.1, sim5.2
)

sim_sce3 <- list(
  sim3.1, sim3.2, sim3.3,
  sim3.1_eta_flat, sim3.2_eta_flat, sim3.3_eta_flat,
  sim3.1_eta_gamma_flat, sim3.2_eta_gamma_flat, sim3.3_eta_gamma_flat,
  sim3.1_all_flat, sim3.2_all_flat, sim3.3_all_flat
)

summary_bin <- get_summary_bin(sim_sces)
summary_bin_sce3 <- get_summary_bin(sim_sce3)

# BCI/coverage lists used by the "plot BCI" section below -- these come out
# of get_summary_bin()'s return value; they used to be referenced (as bare
# `bcis_95`, `coverages_95`, etc.) without ever being extracted here, which
# meant this script could only run at all if a same-named object happened to
# be left over in the environment from a previous session.
bcis_95 <- summary_bin$bcis_95
bcis_90 <- summary_bin$bcis_90
bcis_80 <- summary_bin$bcis_80
bcis_50 <- summary_bin$bcis_50
coverages_95 <- summary_bin$coverages_95
coverages_90 <- summary_bin$coverages_90
coverages_80 <- summary_bin$coverages_80
coverages_50 <- summary_bin$coverages_50

# ------------------------------------------------------------------------------
# Scientific-notation helper for the MSE/WIS summary tables
# ------------------------------------------------------------------------------
# Rather than formatting every cell in scientific notation (e.g. "1.23e-03"),
# each table is shown on a single common power-of-10 scale: all its values
# are multiplied by 10^exponent so they land near order 1, and the exponent
# is reported once in the table title (e.g. "Average MSE (\u00d710\u207b\u00b3)")
# so it can be echoed in the title or an Overleaf figure legend instead of
# cluttering every cell.
.superscript_digits <- c(
  "0" = "\u2070", "1" = "\u00B9", "2" = "\u00B2", "3" = "\u00B3", "4" = "\u2074",
  "5" = "\u2075", "6" = "\u2076", "7" = "\u2077", "8" = "\u2078", "9" = "\u2079"
)

to_superscript <- function(n) {
  neg <- n < 0
  digits <- strsplit(as.character(abs(n)), "")[[1]]
  s <- paste(.superscript_digits[digits], collapse = "")
  if (neg) paste0("\u207B", s) else s
}

# Exponent k such that multiplying x by 10^k brings its median magnitude
# close to order 1 (e.g. x ~ 1e-3 -> k = 3).
compute_common_scale <- function(x) {
  x <- x[is.finite(x) & x != 0]
  if (length(x) == 0) return(0L)
  round(-floor(log10(median(abs(x)))))
}

# Multiply `cols` of `df` by a common power of 10 and return the scaled
# data frame together with a " (\u00d710^-k)" label for the table title --
# the exponent shown is the one needed to recover the *original* values
# from the displayed (scaled) ones.
add_common_scale <- function(df, cols) {
  exponent <- compute_common_scale(unlist(df[cols]))
  df[cols] <- lapply(df[cols], function(col) col * 10^exponent)
  label <- if (exponent == 0) "" else paste0(" (\u00d710", to_superscript(-exponent), ")")
  list(df = df, exponent = exponent, label = label)
}

df_mse_theta <- data.frame(
  # Regime = c(rep("n[0] > n", 12), rep("n[0] < n", 12)),
  Theta0 = c(
    "(0.30, 0.40, 0.50)",
    "(0.20, 0.30, 0.40)",
    "(0.05, 0.15, 0.25)",
    "(0.30, 0.50, 0.40)",
    "(0.20, 0.40, 0.30)",
    "(0.05, 0.25, 0.15)",
    "(0.50, 0.40, 0.30)",
    "(0.40, 0.30, 0.20)",
    "(0.25, 0.15, 0.05)",
    "(0.50, 0.50, 0.50)",
    "(0.30, 0.30, 0.30)",
    "(0.10, 0.10, 0.10)",
    "(0.10, 0.90, 0.50)",
    "(0.50, 0.90, 0.10)"
  ),
  Order       = c(rep("Compatible",    3),
                  rep("Semi-compat.",  3),
                  rep("Incompatible",  3),
                  rep("Neutral",       3),
                  rep("Cancellation effect", 2)
                  ),
  Discrepancy = c(rep(c("None","Small","Large"), 4), rep("Large", 2)),
  NPP         = c(
    lapply(summary_bin$mses, function(x) x[1]) %>% unlist()
  ),
  NPP_SEQ     = c(
    lapply(summary_bin$mses, function(x) x[2]) %>% unlist()
  ),
  ONPP        = c(
    lapply(summary_bin$mses, function(x) x[3]) %>% unlist()
  ),
  ONPP_SEQ    = c(
    lapply(summary_bin$mses, function(x) x[4]) %>% unlist()
  )
)

# df_mse_beta$Regime <- ifelse(df_mse_beta$Regime == "n[0] > n", "n\u2080 > n", "n\u2080 < n")

mse_scale <- add_common_scale(df_mse_theta, c("NPP", "NPP_SEQ", "ONPP", "ONPP_SEQ"))
df_mse_theta <- mse_scale$df

table_mse_theta <- df_mse_theta |>
  gt(groupname_col = "Order") |>
  cols_label(
    Theta0 = "\u03B8\u2080",
    Discrepancy = "Discrepancy",
    NPP         = "NPP",
    NPP_SEQ     = "NPP-SEQ",
    ONPP        = "ONPP",
    ONPP_SEQ    = "ONPP-SEQ"
  ) |>
  fmt_number(
    columns  = c(NPP, NPP_SEQ, ONPP, ONPP_SEQ),
    decimals = 3
  ) |>
  tab_header(
    # title    = "Average MSE by scenario and method",
    # subtitle = "Varying \u03B2",
    title = paste0("Average MSE", mse_scale$label)
  ) |>
  tab_spanner(
    label   = "Method",
    columns = c(NPP, NPP_SEQ, ONPP, ONPP_SEQ)
  ) |>
  tab_style(
    style     = cell_text(weight = "bold"),
    locations = cells_row_groups()
  ) |>
  tab_options(
    row_group.font.weight    = "bold",
    row_group.background.color = "#f0f0f0",
    table.font.size          = 11
  )

# Row-wise gradient as before, plus a bold-white highlight on each row's
# best (smallest) cell -- see highlight_best_per_row() in aux_fun_sim.R.
method_cols <- c("NPP", "NPP_SEQ", "ONPP", "ONPP_SEQ")
table_mse_theta <- highlight_best_per_row(table_mse_theta, df_mse_theta,
                                          columns = method_cols,
                                          palette = viridisLite::mako)

gtsave(table_mse_theta, "results/summary_tables/bin/mse_theta_table_bin.png", zoom = 3)

df_wis_theta <- data.frame(
  # Theta0 = c(
  #   "(0.30, 0.40, 0.50)",
  #   "(0.20, 0.30, 0.40)",
  #   "(0.05, 0.15, 0.25)",
  #   "(0.30, 0.50, 0.40)",
  #   "(0.20, 0.40, 0.30)",
  #   "(0.05, 0.25, 0.15)",
  #   "(0.50, 0.40, 0.30)",
  #   "(0.40, 0.30, 0.20)",
  #   "(0.25, 0.15, 0.05)",
  #   "(0.50, 0.50, 0.50)",
  #   "(0.30, 0.30, 0.30)",
  #   "(0.10, 0.10, 0.10)",
  #   "(0.10, 0.90, 0.50)"
  # ),
  Order       = c(rep("Compatible",    3),
                  rep("Semi-compat.",  3),
                  rep("Incompatible",  3),
                  rep("Neutral",       3),
                  rep("Cancellation effect", 2)
                  ),
  # Discrepancy = c(rep(c("None","Small","Large"), 4), "Large"),
  NPP         = c(
    lapply(summary_bin$wis_list, function(x) x$hat_wis_npp[1]) %>% unlist()
  ),
  NPP_SEQ     = c(
    lapply(summary_bin$wis_list, function(x) x$hat_wis_nppseq[1]) %>% unlist()
  ),
  ONPP        = c(
    lapply(summary_bin$wis_list, function(x) x$hat_wis_onpp[1]) %>% unlist()
  ),
  ONPP_SEQ    = c(
    lapply(summary_bin$wis_list, function(x) x$hat_wis_onppseq[1]) %>% unlist()
  )
)
wis_scale <- add_common_scale(df_wis_theta, c("NPP", "NPP_SEQ", "ONPP", "ONPP_SEQ"))
df_wis_theta <- wis_scale$df

table_wis_theta <- df_wis_theta |>
  gt(groupname_col = "Order") |>
  cols_label(
    # Theta0 = "\u03B8\u2080",
    # Discrepancy = "Discrepancy",
    NPP         = "NPP",
    NPP_SEQ     = "NPP-SEQ",
    ONPP        = "ONPP",
    ONPP_SEQ    = "ONPP-SEQ"
  ) |>
  fmt_number(
    columns  = c(NPP, NPP_SEQ, ONPP, ONPP_SEQ),
    decimals = 3
  ) |>
  tab_header(
    title    = paste0("Average WIS", wis_scale$label)
  ) |>
  tab_spanner(
    label   = "Method",
    columns = c(NPP, NPP_SEQ, ONPP, ONPP_SEQ)
  ) |>
  tab_style(
    style     = cell_text(weight = "bold"),
    locations = cells_row_groups()
  ) |>
  tab_options(
    row_group.font.weight    = "bold",
    row_group.background.color = "#f0f0f0",
    table.font.size          = 11
  )

table_wis_theta <- highlight_best_per_row(table_wis_theta, df_wis_theta,
                                          columns = method_cols,
                                          palette = viridisLite::plasma)
gtsave(table_wis_theta, "results/summary_tables/bin/wis_theta_table_bin.png", zoom = 3)

table_mse_wis_theta <- gt_two_column_layout(list(table_mse_theta, table_wis_theta))

htmltools::save_html(table_mse_wis_theta, "table_mse_wis_theta.html")
webshot2::webshot("table_mse_wis_theta.html",
                  "results/summary_tables/bin/mse_wis_theta_table.png",
                  zoom    = 3,
                  vwidth  = 600,   # reduce until white space disappears
                  vheight = 1200)

# plot BCI
#
# The four blocks below (95/90/80/50%) used to be near-identical
# copy-pasted lapply()s. Each one's "is this the last row" legend check --
# `if (i == length(bcis_95))` -- was copy-pasted along with the rest of the
# block, so the 90/80/50% versions all compared against `length(bcis_95)`
# instead of their own list. That was harmless only because all four lists
# happen to have the same length (14); factoring the block into one helper,
# parametrized by which list it's actually plotting, removes the trap.
build_bci_plot_grid <- function(bci_list, coverage_list, prob) {
  plots <- lapply(seq_along(bci_list), function(i) {
    sce <- as.roman(ceiling(i/3))
    sce <- ifelse(i %% 3 == 1, paste0(sce,".I"), ifelse(i %% 3 == 2, paste0(sce,".II"), paste0(sce,".III")))
    title_left <- ggplot() +
      annotate("text", x = 0, y = 0,
               label = sce,
               angle = 90, size = 8) +
      theme_void()
    plot <- plot_bci(bci_list[[i]], sim_sces[[i]], true_value, prob, coverage_list[[i]], sce, "Bernoulli")
    plot <- (title_left | plot) + plot_layout(widths = c(0.05, 1))
    if (i == length(bci_list)) {
      plot <- plot & theme(legend.position = "bottom")
    }
    return(plot)
  })
  Reduce(`/`, plots)
}

bci_specs <- list(
  list(bcis = bcis_95, coverages = coverages_95, prob = 0.95, file = "results/figures/bin/bci/bci_95_bin.pdf"),
  list(bcis = bcis_90, coverages = coverages_90, prob = 0.90, file = "results/figures/bin/bci/bci_90_bin.pdf"),
  list(bcis = bcis_80, coverages = coverages_80, prob = 0.80, file = "results/figures/bin/bci/bci_80_bin.pdf"),
  list(bcis = bcis_50, coverages = coverages_50, prob = 0.50, file = "results/figures/bin/bci/bci_50_bin.pdf")
)

for (spec in bci_specs) {
  wrap_plots_bci <- build_bci_plot_grid(spec$bcis, spec$coverages, spec$prob)
  ggsave(spec$file, wrap_plots_bci, width = 21, height = 29, dpi = 300, limitsize = F)
}


# ------------------------------------------------------------------------------
# Comparison of scenarios 3 across prior specifications (A)/(B)/(C)/(D)
# ------------------------------------------------------------------------------
#
# sim_sce3 (built above) holds, in order: the three "incompatible order"
# scenarios (3_1/3_2/3_3) under the default prior -- (A) a = 1/2, a' = 1/2,
# alpha = 1/(K+1) -- then the same three scenarios again under three
# progressively flatter priors: (B) eta_flat (a' = 1), (C) eta_gamma_flat
# (a' = 1, alpha = 1), (D) all_flat (a = 1, a' = 1, alpha = 1). So
# summary_bin_sce3$biass/$mses (12 entries, 4 methods each) split into four
# consecutive 3-scenario blocks, one per prior column -- see
# scenarios_bin.R's base_par/flat_variants for where (A)-(D) actually come
# from. This replaces the old flat table (columns labeled
# "Bias_Default"/"MSE_Eta_Flat"/... with no indication of what those prior
# labels actually mean) with the (A)/(B)/(C)/(D) x Bias/MSE two-level
# header and multirow scenario labels used in the paper.
#
# Requires \usepackage{multirow} in the manuscript's LaTeX preamble (for the
# per-scenario multirow label in the first column) -- not something this
# script can add for you.

methods <- c("NPP", "NPP-SEQ", "ONPP", "ONPP-SEQ")

sce3_labels <- data.frame(
  roman       = c("III.I", "III.II", "III.III"),
  theta0      = vapply(sce3_theta0, function(v) paste0("(", paste(formatC(v, format = "f", digits = 2), collapse = ", "), ")"), character(1)),
  discrepancy = c("No discrepancy", "Small discrepancy", "Large discrepancy")
)

# One multirow LaTeX cell per scenario (content only on that block's first
# row, blank on the other 3, exactly as \multirow expects).
scenario_cell <- unlist(lapply(seq_len(nrow(sce3_labels)), function(s) {
  label <- paste0(
    "\\multirow{4}{*}{\\shortstack[l]{Scenario ", sce3_labels$roman[s],
    " \\\\ ", sce3_labels$theta0[s],
    " \\\\ Incompatible order \\\\ ", sce3_labels$discrepancy[s], "}}"
  )
  c(label, "", "", "")
}))

comp_table <- data.frame(
  Theta0  = scenario_cell,
  Prior   = rep(methods, nrow(sce3_labels)),
  Bias_A  = unlist(summary_bin_sce3$biass[1:3]),
  MSE_A   = unlist(summary_bin_sce3$mses[1:3]),
  Bias_B  = unlist(summary_bin_sce3$biass[4:6]),
  MSE_B   = unlist(summary_bin_sce3$mses[4:6]),
  Bias_C  = unlist(summary_bin_sce3$biass[7:9]),
  MSE_C   = unlist(summary_bin_sce3$mses[7:9]),
  Bias_D  = unlist(summary_bin_sce3$biass[10:12]),
  MSE_D   = unlist(summary_bin_sce3$mses[10:12]),
  stringsAsFactors = FALSE
)

comp_table <- comp_table %>%
  mutate(across(c(Bias_A, MSE_A, Bias_B, MSE_B, Bias_C, MSE_C, Bias_D, MSE_D),
                ~ formatC(.x, format = "f", digits = 3)))

# Two-level header: (A)/(B)/(C)/(D) spanners over Bias/MSE sub-columns,
# injected in place of xtable's default single-row column names
# (include.colnames = FALSE below) -- xtable has no native support for a
# grouped/spanning header, so this builds it by hand as raw LaTeX.
addtorow <- list(
  pos = list(-1, -1),
  command = c(
    "\\hline & & \\multicolumn{2}{c}{(A)} & \\multicolumn{2}{c}{(B)} & \\multicolumn{2}{c}{(C)} & \\multicolumn{2}{c}{(D)} \\\\\n",
    "$\\theta_0$ & Prior & Bias & MSE & Bias & MSE & Bias & MSE & Bias & MSE \\\\\n"
  )
)

print(
  xtable::xtable(comp_table, align = c("l", "l", "l", rep("r", 8))),
  add.to.row = addtorow,
  include.colnames = FALSE,
  include.rownames = FALSE,
  hline.after = c(0, 4, 8, 12),
  sanitize.text.function = identity,
  type = "latex"
)

# ------------------------------------------------------------------------------
# Scenario 3 comparison across priors (A)-(D) -- gt-table version, MSE + WIS
# ------------------------------------------------------------------------------
#
# Same visual style as the "Average MSE"/"Average WIS" tables built above
# (row-group header, mako/plasma per-row gradient, bold-white highlight on
# the row's best method -- see highlight_best_per_row() in aux_fun_sim.R),
# but for the prior-sensitivity comparison instead of the discrepancy-order
# comparison: rows are grouped by prior spec -- (A) default, (B) eta_flat,
# (C) eta_gamma_flat, (D) all_flat, see scenarios_bin.R's
# base_par/flat_variants -- and within each group there's one row per
# incompatible-order scenario (III.I/II/III). WIS replaces Bias here: Bias
# isn't on the same scale as MSE/WIS (it can be negative, and its magnitude
# means something different), so it doesn't fit this table's single
# continuous per-row color gradient the way a second error-magnitude metric
# does -- WIS is the same "how good is the interval forecast" companion
# metric used everywhere else in this script.

prior_labels <- c("(A) Default", "(B) η flat", "(C) η, γ flat", "(D) All flat")

sce3_grid <- data.frame(
  Theta0      = rep(sce3_labels$theta0, times = 4),
  Discrepancy = rep(sce3_labels$discrepancy, times = 4),
  Prior       = rep(prior_labels, each = 3)
)

df_mse_sce3 <- cbind(sce3_grid, data.frame(
  NPP      = unlist(lapply(summary_bin_sce3$mses, function(x) x[1])),
  NPP_SEQ  = unlist(lapply(summary_bin_sce3$mses, function(x) x[2])),
  ONPP     = unlist(lapply(summary_bin_sce3$mses, function(x) x[3])),
  ONPP_SEQ = unlist(lapply(summary_bin_sce3$mses, function(x) x[4]))
))

mse_sce3_scale <- add_common_scale(df_mse_sce3, c("NPP", "NPP_SEQ", "ONPP", "ONPP_SEQ"))
df_mse_sce3 <- mse_sce3_scale$df

table_mse_sce3 <- df_mse_sce3 |>
  gt(groupname_col = "Prior") |>
  cols_label(
    Theta0      = "θ₀",
    Discrepancy = "Discrepancy",
    NPP         = "NPP",
    NPP_SEQ     = "NPP-SEQ",
    ONPP        = "ONPP",
    ONPP_SEQ    = "ONPP-SEQ"
  ) |>
  fmt_number(
    columns  = c(NPP, NPP_SEQ, ONPP, ONPP_SEQ),
    decimals = 3
  ) |>
  tab_header(
    title = paste0("Average MSE", mse_sce3_scale$label)
  ) |>
  tab_spanner(
    label   = "Method",
    columns = c(NPP, NPP_SEQ, ONPP, ONPP_SEQ)
  ) |>
  tab_style(
    style     = cell_text(weight = "bold"),
    locations = cells_row_groups()
  ) |>
  tab_options(
    row_group.font.weight    = "bold",
    row_group.background.color = "#f0f0f0",
    table.font.size          = 11
  )

table_mse_sce3 <- highlight_best_per_row(table_mse_sce3, df_mse_sce3,
                                         columns = method_cols,
                                         palette = viridisLite::mako)

gtsave(table_mse_sce3, "results/summary_tables/bin/mse_sce3_table_bin.png", zoom = 3)

df_wis_sce3 <- cbind(sce3_grid, data.frame(
  NPP      = unlist(lapply(summary_bin_sce3$wis_list, function(x) x$hat_wis_npp[1])),
  NPP_SEQ  = unlist(lapply(summary_bin_sce3$wis_list, function(x) x$hat_wis_nppseq[1])),
  ONPP     = unlist(lapply(summary_bin_sce3$wis_list, function(x) x$hat_wis_onpp[1])),
  ONPP_SEQ = unlist(lapply(summary_bin_sce3$wis_list, function(x) x$hat_wis_onppseq[1]))
))

wis_sce3_scale <- add_common_scale(df_wis_sce3, c("NPP", "NPP_SEQ", "ONPP", "ONPP_SEQ"))
df_wis_sce3 <- wis_sce3_scale$df
df_wis_sce3 <- df_wis_sce3 %>%
  select(Prior, NPP, NPP_SEQ, ONPP, ONPP_SEQ)

table_wis_sce3 <- df_wis_sce3 |>
  gt(groupname_col = "Prior") |>
  cols_label(
    # Theta0      = "θ₀",
    # Discrepancy = "Discrepancy",
    NPP         = "NPP",
    NPP_SEQ     = "NPP-SEQ",
    ONPP        = "ONPP",
    ONPP_SEQ    = "ONPP-SEQ"
  ) |>
  fmt_number(
    columns  = c(NPP, NPP_SEQ, ONPP, ONPP_SEQ),
    decimals = 3
  ) |>
  tab_header(
    title = paste0("Average WIS", wis_sce3_scale$label)
  ) |>
  tab_spanner(
    label   = "Method",
    columns = c(NPP, NPP_SEQ, ONPP, ONPP_SEQ)
  ) |>
  tab_style(
    style     = cell_text(weight = "bold"),
    locations = cells_row_groups()
  ) |>
  tab_options(
    row_group.font.weight    = "bold",
    row_group.background.color = "#f0f0f0",
    table.font.size          = 11
  )

table_wis_sce3 <- highlight_best_per_row(table_wis_sce3, df_wis_sce3,
                                         columns = method_cols,
                                         palette = viridisLite::plasma)

gtsave(table_wis_sce3, "results/summary_tables/bin/wis_sce3_table_bin.png", zoom = 3)

table_mse_wis_sce3 <- gt_two_column_layout(list(table_mse_sce3, table_wis_sce3))

htmltools::save_html(table_mse_wis_sce3, "table_mse_wis_sce3.html")
webshot2::webshot("table_mse_wis_sce3.html",
                  "results/summary_tables/bin/mse_wis_sce3_table.png",
                  zoom    = 3,
                  vwidth  = 630
                )

# ------------------------------------------------------------------------------
# Plot boxplots
# ------------------------------------------------------------------------------

plots_bin <- lapply(1:length(sim_sces), plot_sce_bin)
lapply(1:4, function(i) {
  sce_plots <- lapply(1:3, function(j) {
    plots_bin[[(i-1)*3 + j]]
  })
  combined_plot <- Reduce(`/`, sce_plots)
  sce <- as.roman(i)
  ggsave(paste0("results/figures/bin/boxplot_sce_", sce, "_bin.pdf"), combined_plot, width = 10, height = 12, dpi = 300)
})


# ------------------------------------------------------------------------------
# Count how many times NPP and NPP-SEQ get the order right
# ------------------------------------------------------------------------------

order_list <- list(c(1,3),
                   c(1,3),
                   c(1,3),
                   c(1,2),
                   c(1,2),
                   c(1,2),
                   c(3,1),
                   c(3,1),
                   c(3,1)
)

order_counts <- lapply(1:length(order_list), function(i) {
  sim <- sim_sces[[i]]
  min_max <- order_list[[i]]
  delta_npp <- lapply(sim$delta, function(x) x$delta_npp)
  delta_nppseq <- lapply(sim$delta, function(x) x$delta_nppseq)
  count_npp <- lapply(delta_npp, function(x) {
    (apply(x, 1, which.min) == min_max[1]) &
      (apply(x, 1, which.max) == min_max[2])
  })
  avg_count_npp <- unlist(count_npp) %>% mean()
  count_nppseq <- lapply(delta_nppseq, function(x) {
    (apply(x, 1, which.min) == min_max[1]) &
      (apply(x, 1, which.max) == min_max[2])
  })
  avg_count_nppseq <- unlist(count_nppseq) %>% mean()
  return(list(
    npp = avg_count_npp,
    nppseq = avg_count_nppseq
  ))
})
order_counts_df <- data.frame(
  scenario = unlist(lapply(1:length(order_list), 
                           function(i) {
                             sce <- as.roman(ceiling(i/3))
                             sce <- ifelse(i %% 3 == 1, paste0(sce,".I"), ifelse(i %% 3 == 2, paste0(sce,".II"), paste0(sce,".III")))
                             rep(sce,2)
                           })
  ),
  model = rep(c("NPP", "NPP-SEQ"), length(order_list)),
  order_count = unlist(lapply(order_counts, function(x) c(x$npp, x$nppseq)))
)    
print(xtable::xtable(order_counts_df), 
      type = "latex", include.rownames = FALSE)

# ------------------------------------------------------------------------------
# Compute how badly ONPP and ONPP-SEQ estimate the order
# ------------------------------------------------------------------------------

order_error_list <- lapply(1:12, function(i) {
  sim <- sim_sces[[i]]
  delta_onpp <- lapply(sim$delta, function(x) x$delta_onpp)
  delta_onppseq <- lapply(sim$delta, function(x) x$delta_onppseq)
  avg_dist1_onpp <- unlist(lapply(delta_onpp, function(x) {
    apply(x, 1, function(y) {
      y[2] - y[1]
    })
  })) %>% mean()
  avg_dist2_onpp <- unlist(lapply(delta_onpp, function(x) {
    apply(x, 1, function(y) {
      y[3] - y[2]
    })
  })) %>% mean()
  avg_dist1_onppseq <- unlist(lapply(delta_onppseq, function(x) {
    apply(x, 1, function(y) {
      y[2] - y[1]
    })
  })) %>% mean()
  avg_dist2_onppseq <- unlist(lapply(delta_onppseq, function(x) {
    apply(x, 1, function(y) {
      y[3] - y[2]
    })
  })) %>% mean()
  return(list(
    onpp = c(avg_dist1_onpp, avg_dist2_onpp),
    onppseq = c(avg_dist1_onppseq, avg_dist2_onppseq)
  ))
})

order_error_df <- data.frame(
  scenario = unlist(lapply(1:12, 
                           function(i) {
                             sce <- as.roman(ceiling(i/3))
                             sce <- ifelse(i %% 3 == 1, paste0(sce,".I"), ifelse(i %% 3 == 2, paste0(sce,".II"), paste0(sce,".III")))
                             rep(sce,2)
                           })
  ),
  model = rep(c("ONPP", "ONPP-SEQ"), 12),
  dist_med_min = unlist(lapply(order_error_list, function(x) c(x$onpp[1], x$onppseq[1]))),
  dist_max_med = unlist(lapply(order_error_list, function(x) c(x$onpp[2], x$onppseq[2])))
)  
print(xtable::xtable(order_error_df %>%
                       mutate(
                         dist_med_min = formatC(dist_med_min, format = "f", digits = 3),
                         dist_max_med = formatC(dist_max_med, format = "f", digits = 3)
                       )),
      type = "latex", include.rownames = FALSE)


# ------------------------------------------------------------------------------
# Plot distributions for a given simulation
# ------------------------------------------------------------------------------

sim <- 5  # choose which simulation to plot

sce3_priors <- list(
  list(
    sim3.1, sim3.1_eta_flat, sim3.1_eta_gamma_flat, sim3.1_all_flat
  ),
  list(
    sim3.2, sim3.2_eta_flat, sim3.2_eta_gamma_flat, sim3.2_all_flat
  ),
  list(
     sim3.3, sim3.3_eta_flat, sim3.3_eta_gamma_flat, sim3.3_all_flat
  )
)

# plot_sce3_eta()/plot_sce3_theta() were rewritten to use facet_grid()/
# facet_wrap() internally (see aux_fun_sim.R) -- each call now returns one
# complete ggplot per scenario (with its own (A)-(D) strip labels and a
# single legend), not a list of per-panel plots, so all that's left to do
# here is stack the 3 scenarios and collect their (now-identical) legends
# into one shared one at the bottom instead of repeating it 3 times.
plots_sce3_eta <- lapply(1:3, function(i) {
  plot_sce3_eta(i, sim = sim, sim_sces = sce3_priors)
})
ggsave("results/figures/bin/post_eta_sce_III_priors_bin.pdf",
       wrap_plots(plots_sce3_eta, ncol = 1) +
         plot_layout(guides = "collect") &
         theme(legend.position = "bottom"),
       width = 10, height = 13, dpi = 300,
       device = cairo_pdf
      )

plots_sce3_theta <- lapply(1:3, function(i) {
  plot_sce3_theta(i, sim = sim, sim_sces = sce3_priors)
})
ggsave("results/figures/bin/post_theta_sce_III_priors_bin.pdf",
       wrap_plots(plots_sce3_theta, ncol = 1) +
         plot_layout(guides = "collect") &
         theme(legend.position = "bottom"),
       width = 10, height = 8, dpi = 300,
       device = cairo_pdf  
      )


# NOTE: this originally tried to compare prior-only draws across all four
# scenario III prior specs (default / eta_flat / eta_gamma_flat / all_flat),
# via `samples_prior_sce3.1`, `samples_prior_sce3.1_eta_flat` and
# `samples_prior_sce3.1_eta_gamma_flat`. Those three were never generated by
# simulations_bin.R (only the "all flat" prior-only draws are ever sampled --
# see samples_prior_sce3.1_all_flat above), so referencing them here always
# errored. Per decision, this now plots only the one prior spec that's
# actually available; extend simulations_bin.R to sample the other three
# scenario III prior specs (prior-only, post = 0) if the full 4-way
# comparison is needed later.
samples_sce3_eta_prior <- list(
    samples_prior_sce3.1_all_flat
)

plots_eta_prior <- plot_eta_prior(samples_sce3_eta_prior)
ggsave("results/figures/bin/prior_eta_priors_bin.pdf",
       plots_eta_prior, width = 10, height = 6, dpi = 300)

samples_sce3_theta_prior <- list(
    rbeta(8000, scenarios[["3_1_all_flat"]]$a, scenarios[["3_1_all_flat"]]$b)
)
plots_theta_prior <- plot_theta_prior(samples_sce3_theta_prior)
plots_theta_prior
ggsave("results/figures/bin/prior_theta_priors_bin.pdf",
       plots_theta_prior, width = 10, height = 3, dpi = 300)


#------------------------------------------------------------------------------
# Cancellation effect scenario
#------------------------------------------------------------------------------

plots_eta_post_canc_eff1_default <- plot_eta_post(sim_sces[[13]], 1)
plots_eta_post_canc_eff1_seq <- plot_eta_post(sim_sces[[13]], 1, seq_norm = T)
plots_eta_post_canc_eff1 <- plots_eta_post_canc_eff1_default / plots_eta_post_canc_eff1_seq
ggsave("results/figures/bin/post_eta_canc_eff1_bin.png", 
       plots_eta_post_canc_eff1, width = 10, height = 6, dpi = 300)
plots_eta_post_canc_eff2_default <- plot_eta_post(sim_sces[[14]], 1)
plots_eta_post_canc_eff2_seq <- plot_eta_post(sim_sces[[14]], 1, seq_norm = T)
plots_eta_post_canc_eff2 <- plots_eta_post_canc_eff2_default / plots_eta_post_canc_eff2_seq
ggsave("results/figures/bin/post_eta_canc_eff2_bin.png",
       plots_eta_post_canc_eff2, width = 10, height = 6, dpi = 300)
