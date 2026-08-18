library(parallel)
library(dplyr)
library(tidyverse)

#-------------------------------------------------------------------------------
# Utility functions
#-------------------------------------------------------------------------------

# function to flatten the data
flatten_data <- function(data) {
  flat_data <- lapply(data, as.vector)
  dims_data <- t(sapply(data, dim))
  data_flat <- unlist(flat_data)
  starting_index <- cumsum(c(1, head(dims_data[,1] * dims_data[,2], -1)))
  return(list(data_flat = data_flat, dims_data = dims_data, starting_index = starting_index))
}


# function to prepare data for Stan
prepare_data_stan <- function(data) {
  y0 <- data %>%
    filter(grepl("hist", data_id)) %>%
    group_split(replicate) %>%
    lapply(function(df) df %>% group_split(data_id) %>%
              lapply(function(df) as.matrix(df$y))
          )
  y <- data %>%
    filter(grepl("curr", data_id)) %>%
    group_split(replicate) %>%
    lapply(function(df) as.matrix(df$y))
  X0 <- data %>%
    filter(grepl("hist", data_id)) %>%
    group_split(replicate) %>%
    lapply(function(df) df %>% group_split(data_id) %>%
              lapply(function(df) cbind(1, as.matrix(df %>% select(starts_with("X")))))
          )
  X <- data %>%
    filter(grepl("curr", data_id)) %>%
    group_split(replicate) %>%
    lapply(function(df)
            cbind(1, as.matrix(df %>% select(starts_with("X"))))
          )

  # Flatten the data
  flattened_X0 <- lapply(X0, function(x) flatten_data(x))
  flattened_y0 <- lapply(y0, function(x) flatten_data(x))
  out <- lapply(1:length(flattened_X0), function(i) {
    list(X0_flat = flattened_X0[[i]]$data_flat,
         dims_X0 = flattened_X0[[i]]$dims_data,
         startid_X0 = flattened_X0[[i]]$starting_index,
         y0_flat = flattened_y0[[i]]$data_flat,
         dims_y0 = flattened_y0[[i]]$dims_data,
         startid_y0 = flattened_y0[[i]]$starting_index,
         X = X[[i]],
         y = y[[i]]
        )
  })
  return(out)
}

#-------------------------------------------------------------------------------
# Scenario-building helpers
#-------------------------------------------------------------------------------
#
# This file used to be ~1080 lines: the same "read 3 CSVs, wrap each in a
# scenario list with an identical `par`, swap two historical-arm labels to
# build the reordered variants, qs_save() 12-14 files" recipe copy-pasted
# once per (varying beta/sigma) x (n0 >= n / n0 <= n) combination. Besides
# the sheer duplication, that copy-pasting let a real bug slip in: the
# "Order II" block's swap always checked for `data_id` values "hist2"/
# "hist3" (no underscore), while the actual column values are "hist_1"/
# "hist_2"/"hist_3" (confirmed by the very next block, "Order III", which
# correctly uses the underscored names). Because "hist2"/"hist3" never
# matched anything, `idx2`/`idx3` were always all-FALSE and the relabeling
# was a silent no-op -- every "Order II" (sceII_*) scenario was actually
# built from unmodified Order I data, in all four sections of the original
# script. If you've already run this script and used its sceII_* outputs
# (directly, or via simulations_lm.r/analysis results downstream), those
# results reflect Order I data, not the intended reordered arms, and should
# be regenerated.

lm_default_par <- list(p = 3, a = 2, b = 1, V0 = diag(3), mu0 = rep(0, 3),
                        al = 1, bl = 1, alpha = rep(1, 4))

# Swaps two historical-arm data_id labels in place (e.g. "hist_2"/"hist_3"),
# used to build the "reordered arms" scenario variants from the baseline
# data without needing separate CSV files for every reordering.
swap_data_id_labels <- function(df, label_a, label_b) {
  idx_a <- df$data_id == label_a
  idx_b <- df$data_id == label_b
  df$data_id[idx_a] <- label_b
  df$data_id[idx_b] <- label_a
  df
}

make_scenario <- function(data, par = lm_default_par) {
  list(par = par, data = prepare_data_stan(data))
}

# Builds and saves the 12 sce<Order>_<Congruence>_<suffix>.qs2 files (4
# orders x 3 congruence levels) for one varying-type/n0-relation
# combination (e.g. csv_type = "beta", n0_rel = "ge" covers the
# "varying_beta_n0_ge_n" section). Order I is the data as-is; Order II
# swaps the hist_2/hist_3 labels; Order III swaps hist_1/hist_3; Order IV
# uses the separate "_neutral_" CSVs rather than a label swap. Congruence
# I/II/III = high/small/no, per the "_high_cong_"/"_small_cong_"/
# "_no_cong_" CSV naming.
build_lm_order_scenarios <- function(csv_type, n0_rel, num_sim, par = lm_default_par) {
  suffix <- sprintf("varying_%s_n0_%s_n", csv_type, n0_rel)
  congruence <- c(I = "high_cong", II = "small_cong", III = "no_cong")

  read_one <- function(cong_token, neutral = FALSE) {
    neutral_part <- if (neutral) "neutral_" else ""
    path <- sprintf("../../simulated-data-regression/data/sim_data_%s_%svarying_%s_p3_n0_%s_n.csv",
                     cong_token, neutral_part, csv_type, n0_rel)
    read_csv(path) %>% filter(replicate <= num_sim)
  }

  order_I   <- lapply(congruence, read_one)
  order_II  <- lapply(order_I, swap_data_id_labels, label_a = "hist2", label_b = "hist3")
  order_III <- lapply(order_I, swap_data_id_labels, label_a = "hist1", label_b = "hist3")
  order_IV  <- lapply(congruence, read_one, neutral = TRUE)

  orders <- list(I = order_I, II = order_II, III = order_III, IV = order_IV)

  scenarios <- list()
  for (order_key in names(orders)) {
    for (cong_key in names(congruence)) {
      key <- paste0(order_key, "_", cong_key)
      scenarios[[key]] <- make_scenario(orders[[order_key]][[cong_key]], par)
      file <- sprintf("results/sim_data/lm/sce%s_%s.qs2", key, suffix)
      message("prepare_data_lm: ", key, " (", suffix, ") -> ", file)
      qs2::qs_save(scenarios[[key]], file = file)
    }
  }
  invisible(scenarios)
}

#-------------------------------------------------------------------------------
# Generate scenarios
#-------------------------------------------------------------------------------

num_sim <- 200

## varying beta, n0 >= n -- small size (Order I only: high vs no congruence;
## no small-congruence, reordered, or neutral variant exists at this size)
data_sce1.1_small <- read_csv("../../simulated-data-regression/data/sim_data_high_cong_varying_beta_p3_n0_ge_n_small_size.csv") %>%
  filter(replicate <= num_sim)
data_sce1.3_small <- read_csv("../../simulated-data-regression/data/sim_data_incong_varying_beta_p3_n0_ge_n_small_size.csv") %>%
  filter(replicate <= num_sim)

qs2::qs_save(make_scenario(data_sce1.1_small),
  file = "results/sim_data/lm/sceI_I_varying_beta_n0_ge_n_small_size.qs2"
)
qs2::qs_save(make_scenario(data_sce1.3_small),
  file = "results/sim_data/lm/sceI_III_varying_beta_n0_ge_n_small_size.qs2"
)

## varying beta, n0 >= n
build_lm_order_scenarios("beta", "ge", num_sim)

## varying beta, n0 < n
build_lm_order_scenarios("beta", "le", num_sim)

## varying sigma, n0 > n
build_lm_order_scenarios("sigma", "ge", num_sim)

## varying sigma, n0 < n
build_lm_order_scenarios("sigma", "le", num_sim)
