# =============================================================================
# PLOTTING step for the NPP/ONPP/HLNNPP real-data comparison (lupus, HOVON
# -- see real_data_sample.R's header for what each dataset is and why the
# third, NSCLC, example from Script-Example.R isn't included). Builds every
# figure/table from the cache real_data_sample.R writes -- NO Stan sampling
# happens here. Run real_data_sample.R first if that cache doesn't exist.
#
# "No borrowing" below always means the reference posterior real_data_
# sample.R cached (Beta fit to the FULL trial's own sample size, per
# Script-Example.R's lupus_t/HOVON_t) -- i.e. "what you'd conclude without
# any historical data at all, if you'd waited for the full trial."
#
# Run from the project root.
# =============================================================================

library(dplyr)
library(tidyverse)
library(patchwork)
library(gt)
library(qs2)

source("code/aux_fun_sim.R")

cache_file <- "results/samples/real_data/real_data_fits.qs2"
if (!file.exists(cache_file)) {
  stop("Cache not found at ", cache_file, " -- run code/real_data/",
       "real_data_sample.R first (it produces this cache).")
}
cached   <- qs2::qs_read(cache_file)
datasets <- cached$datasets
fits     <- cached$fits

dataset_labels <- setNames(vapply(datasets, `[[`, character(1), "label"), names(datasets))
dataset_order  <- names(datasets)

method_labels <- c(npp = "NPP", onpp = "ONPP", logit = "HLNNPP")
# Same palette as code/poster/prior_comparison_mse_ratio_forest_plot.R
# (ONPP = "#0072B2", HLN-NPP = "#009E73"). NPP itself is never given a
# point color there -- it's the baseline, drawn as a plain black
# geom_vline(xintercept = 1) -- so it gets plain black here too.
method_colors <- c("NPP" = "black", "ONPP" = "#0072B2", "HLNNPP" = "#009E73")
method_order  <- c("npp", "onpp", "logit")

# Empirical HPD interval -- same definition as Script-Example.R's emp.hpd().
emp_hpd <- function(x, conf = 0.95) {
  conf <- min(conf, 1 - conf)
  n  <- length(x)
  nn <- round(n * conf)
  x  <- sort(x)
  xx <- x[(n - nn + 1):n] - x[1:nn]
  nnn <- which(xx == min(xx))[1]
  c(lower = x[nnn], upper = x[n - nn + nnn])
}

#-------------------------------------------------------------------------------
# theta draws, 3 priors + "no borrowing" reference, long format (the
# combined density+forest figure below is built from this).
#-------------------------------------------------------------------------------

theta_long <- bind_rows(lapply(dataset_order, function(key) {
  fit <- fits[[key]]
  bind_rows(
    bind_rows(lapply(method_order, function(mkey) {
      data.frame(value = fit[[paste0("theta_", mkey)]], method = method_labels[[mkey]])
    })),
    data.frame(value = fit$theta_ref, method = "No borrowing")
  ) %>% mutate(dataset = dataset_labels[[key]])
}))
theta_long$method  <- factor(theta_long$method, levels = c(method_labels[method_order], "No borrowing"))
theta_long$dataset <- factor(theta_long$dataset, levels = unname(dataset_labels[dataset_order]))

theta_colors <- c(method_colors, "No borrowing" = "grey40")
theta_linetypes <- c(setNames(rep("solid", 3), method_labels[method_order]), "No borrowing" = "dashed")

#-------------------------------------------------------------------------------
# eta (per-historical-dataset discount) boxplots -- one panel per historical
# dataset k = 1..K, 3 methods side by side, per dataset (lupus: K = 2,
# HOVON: K = 4) -- the real-data analogue of Script-Example.R's box1.1/
# box1.2/box2.1-box2.4, but with HLNNPP added as a third method (that
# script only had the 2 NPP/2 ONPP variants).
#-------------------------------------------------------------------------------

make_eta_boxplot_panel <- function(fit, k, dataset_key) {
  df <- bind_rows(lapply(method_order, function(mkey) {
    mat <- fit[[paste0("eta_", mkey)]]
    data.frame(value = as.numeric(mat[, k]), method = method_labels[[mkey]])
  }))
  df$method <- factor(df$method, levels = method_labels[method_order])

  ggplot(df, aes(x = method, y = value, fill = method)) +
    theme_bw() +
    stat_boxplot(geom = "errorbar", coef = 0.5, width = 0.25) +
    geom_boxplot(outlier.colour = "black", outlier.shape = NA,
                 width = 0.4, outlier.size = 0.3, coef = 0.5) +
    stat_summary(fun = mean, colour = "yellow3", geom = "point",
                  shape = 18, size = 3, show.legend = FALSE) +
    scale_fill_manual(values = method_colors, guide = "none") +
    labs(x = NULL, y = bquote(eta[.(k)])) +
    theme(axis.title.y = element_text(angle = -270, face = "bold", size = 12, vjust = 0.5))
}

for (key in dataset_order) {
  fit <- fits[[key]]
  K   <- datasets[[key]]$K
  panels <- lapply(seq_len(K), function(k) make_eta_boxplot_panel(fit, k, key))
  p_eta <- wrap_plots(panels, nrow = 1)
  ggsave(sprintf("results/figures/real_data/%s_eta_boxplot.pdf", key), p_eta,
         width = 2.6 * K, height = 4, device = cairo_pdf)
}

#-------------------------------------------------------------------------------
# theta posterior density + a forest-style mean/95% HPD summary, combined
# side by side into ONE figure per dataset (3 priors + "no borrowing"
# reference in both panels).
#-------------------------------------------------------------------------------

forest_df <- theta_long %>%
  group_by(dataset, method) %>%
  summarise(mean = mean(value),
            lower = emp_hpd(value)[["lower"]],
            upper = emp_hpd(value)[["upper"]],
            .groups = "drop")
forest_df$method <- factor(forest_df$method, levels = rev(c(method_labels[method_order], "No borrowing")))

for (key in dataset_order) {
  df_density <- theta_long %>% filter(dataset == dataset_labels[[key]])
  df_forest  <- forest_df  %>% filter(dataset == dataset_labels[[key]])

  p_theta_density <- ggplot(df_density, aes(x = value, color = method, linetype = method)) +
    geom_density(linewidth = 0.9) +
    scale_color_manual(name = NULL, values = theta_colors) +
    scale_linetype_manual(name = NULL, values = theta_linetypes) +
    labs(x = expression(theta), y = NULL) +
    .sce3_theme() +
    theme(legend.position = "bottom")

  p_forest <- ggplot(df_forest, aes(x = mean, y = method, xmin = lower, xmax = upper, color = method)) +
    geom_pointrange(size = 0.7, fatten = 3) +
    scale_color_manual(values = theta_colors, guide = "none") +
    labs(x = expression(theta ~ "(posterior mean, 95% HPD)"), y = NULL) +
    .sce3_theme()

  p_theta_combined <- p_theta_density + p_forest

  ggsave(sprintf("results/figures/real_data/%s_theta_density_forest.pdf", key), p_theta_combined,
         width = 11, height = 4.5, device = cairo_pdf)
}

#-------------------------------------------------------------------------------
# HLNNPP hyperparameters (mu, tau, rho) -- diagnostic, same spirit as
# logitnorm_aux_plot.R's hyperparameter densities, faceted by dataset
# instead of by scenario.
#-------------------------------------------------------------------------------

hyper_labels <- c(mu = "μ", tau = "τ", rho = "ρ")

for (key in dataset_order) {
  mat <- fits[[key]]$aux_logit
  df_ds <- bind_rows(lapply(c("mu", "tau", "rho"), function(pname) {
    data.frame(value = as.numeric(mat[, pname]), hyperparam = hyper_labels[[pname]])
  }))
  df_ds$hyperparam <- factor(df_ds$hyperparam, levels = unname(hyper_labels))

  p_hyper <- ggplot(df_ds, aes(x = value, color = hyperparam)) +
    geom_density(linewidth = 0.9) +
    facet_wrap(~ hyperparam, nrow = 1, scales = "free", labeller = label_parsed) +
    scale_color_manual(values = c("#2C5F8A", "#D55E00", "#009E73"), guide = "none") +
    labs(x = NULL, y = NULL) +
    .sce3_theme()

  ggsave(sprintf("results/figures/real_data/%s_hlnnpp_hyperparameters.pdf", key), p_hyper,
         width = 8, height = 3.5, device = cairo_pdf)
}

#-------------------------------------------------------------------------------
# Summary table: theta posterior mean/median/95% HPD by method and dataset.
#-------------------------------------------------------------------------------

summary_df <- theta_long %>%
  group_by(dataset, method) %>%
  summarise(Mean = mean(value), Median = median(value),
            HPD_low = emp_hpd(value)[["lower"]], HPD_high = emp_hpd(value)[["upper"]],
            .groups = "drop")

table_theta_summary <- summary_df |>
  gt(groupname_col = "dataset") |>
  cols_label(method = "Prior") |>
  fmt_number(columns = c(Mean, Median, HPD_low, HPD_high), decimals = 3) |>
  cols_merge(columns = c(HPD_low, HPD_high), pattern = "({1}, {2})") |>
  cols_label(HPD_low = "95% HPD") |>
  tab_header(title = "Posterior summary of the current-study response rate",
             subtitle = "By prior and dataset") |>
  tab_style(style = cell_text(weight = "bold"), locations = cells_row_groups()) |>
  tab_options(row_group.font.weight = "bold", row_group.background.color = "#f0f0f0",
              table.font.size = 12)

gtsave(table_theta_summary, "results/summary_tables/real_data/theta_summary_table.png", zoom = 3)

message("Done. Figures written to results/figures/real_data/, table to ",
        "results/summary_tables/real_data/theta_summary_table.png")
