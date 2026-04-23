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
# Generate scenarios
#-------------------------------------------------------------------------------

num_sim <- 200
data_high_cong <- read_csv("../simulated-data-regression/data/sim_data_high_cong_p3.csv") %>%
  filter(replicate <= num_sim)
data_small_cong <- read_csv("../simulated-data-regression/data/sim_data_small_cong_p3.csv") %>%
  filter(replicate <= num_sim)
data_no_cong <- read_csv("../simulated-data-regression/data/sim_data_no_cong_p3.csv") %>%
  filter(replicate <= num_sim)
data_list <- list(
  data_high_cong,
  data_small_cong,
  data_no_cong
)


sce1.1 <- list(
  par = list(
    p = 3,
    a = 2,
    b = 1,
    V0 = diag(3),
    mu0 = rep(0, 3),
    tilde_a = 1,
    tilde_b = 1,
    alpha = rep(1, 4)
  ),
  data = prepare_data_stan(data_list[[1]])
)
sce1.2 <- list(
  par = list(
    p = 3,
    a = 2,
    b = 1,
    V0 = diag(3),
    mu0 = rep(0, 3),
    tilde_a = 1,
    tilde_b = 1,
    alpha = rep(1, 4)
  ),
  data = prepare_data_stan(data_list[[2]])
)
sce1.3 <- list(
  par = list(
    p = 3,
    a = 2,
    b = 1,
    V0 = diag(3),
    mu0 = rep(0, 3),
    tilde_a = 1,
    tilde_b = 1,
    alpha = rep(1, 4)
  ),
  data = prepare_data_stan(data_list[[3]])
)

# change order
data_sce2.1 <- lapply(data_list, function(df) { 
  df$data_id[df$data_id == "hist_3"] <- "temp"
  df$data_id[df$data_id == "hist_2"] <- "hist_3"
  df$data_id[df$data_id == "temp"] <- "hist_2"
  return(df)
})

sce2.1 <- list(
  par = list(
    p = 3,
    a = 2,
    b = 1,
    V0 = diag(3),
    mu0 = rep(0, 3),
    tilde_a = 1,
    tilde_b = 1,
    alpha = rep(1, 4)
  ),
  data = prepare_data_stan(data_sce2.1[[1]])
)
sce2.2 <- list(
  par = list(
    p = 3,
    a = 2,
    b = 1,
    V0 = diag(3),
    mu0 = rep(0, 3),
    tilde_a = 1,
    tilde_b = 1,
    alpha = rep(1, 4)
  ),
  data = prepare_data_stan(data_sce2.1[[2]])
)
sce2.3 <- list(
  par = list(
    p = 3,
    a = 2,
    b = 1,
    V0 = diag(3),
    mu0 = rep(0, 3),
    tilde_a = 1,
    tilde_b = 1,
    alpha = rep(1, 4)
  ),
  data = prepare_data_stan(data_sce2.1[[3]])
)

data_sce3.1 <- lapply(data_list, function(df) { 
  df$data_id[df$data_id == "hist_3"] <- "temp"
  df$data_id[df$data_id == "hist_1"] <- "hist_3"
  df$data_id[df$data_id == "temp"] <- "hist_1"
  return(df)
})
sce3.1 <- list(
  par = list(
    p = 3,
    a = 2,
    b = 1,
    V0 = diag(3),
    mu0 = rep(0, 3),
    tilde_a = 1,
    tilde_b = 1,
    alpha = rep(1, 4)
  ),
  data = prepare_data_stan(data_sce3.1[[1]])
)
sce3.2 <- list(
  par = list(
    p = 3,
    a = 2,
    b = 1,
    V0 = diag(3),
    mu0 = rep(0, 3),
    tilde_a = 1,
    tilde_b = 1,
    alpha = rep(1, 4)
  ),
  data = prepare_data_stan(data_sce3.1[[2]])
)
sce3.3 <- list(
  par = list(
    p = 3,
    a = 2,
    b = 1,
    V0 = diag(3),
    mu0 = rep(0, 3),
    tilde_a = 1,
    tilde_b = 1,
    alpha = rep(1, 4)
  ),
  data = prepare_data_stan(data_sce3.1[[3]])
)

data_high_cong_neutral <- read_csv("../simulated-data-regression/data/sim_data_high_cong_neutral_p3.csv") %>%
  filter(replicate <= num_sim)
data_small_cong_neutral <- read_csv("../simulated-data-regression/data/sim_data_small_cong_neutral_p3.csv") %>%
  filter(replicate <= num_sim)
data_no_cong_neutral <- read_csv("../simulated-data-regression/data/sim_data_no_cong_neutral_p3.csv") %>%
  filter(replicate <= num_sim)

sce4.1 <- list(
  par = list(
    p = 3,
    a = 2,
    b = 1,
    V0 = diag(3),
    mu0 = rep(0, 3),
    tilde_a = 1,
    tilde_b = 1,
    alpha = rep(1, 4)
  ),
  data = prepare_data_stan(data_high_cong_neutral)
)
sce4.2 <- list(
  par = list(
    p = 3,
    a = 2,
    b = 1,
    V0 = diag(3),
    mu0 = rep(0, 3),
    tilde_a = 1,
    tilde_b = 1,
    alpha = rep(1, 4)
  ),
  data = prepare_data_stan(data_small_cong_neutral)
)
sce4.3 <- list(
  par = list(
    p = 3,
    a = 2,
    b = 1,
    V0 = diag(3),
    mu0 = rep(0, 3),
    tilde_a = 1,
    tilde_b = 1,
    alpha = rep(1, 4)
  ),
  data = prepare_data_stan(data_no_cong_neutral)
)

# save scenarios
qs2::qs_save(sce1.1,
  file = "results/sim_data/lm/sceI_I.qs2"
)
qs2::qs_save(sce1.2,
  file = "results/sim_data/lm/sceI_II.qs2"
)
qs2::qs_save(sce1.3,
  file = "results/sim_data/lm/sceI_III.qs2"
)
qs2::qs_save(sce2.1,
  file = "results/sim_data/lm/sceII_I.qs2"
)
qs2::qs_save(sce2.2,
  file = "results/sim_data/lm/sceII_II.qs2"
)
qs2::qs_save(sce2.3,
  file = "results/sim_data/lm/sceII_III.qs2"
)
qs2::qs_save(sce3.1,
  file = "results/sim_data/lm/sceIII_I.qs2"
)
qs2::qs_save(sce3.2,
  file = "results/sim_data/lm/sceIII_II.qs2"
)
qs2::qs_save(sce3.3,
  file = "results/sim_data/lm/sceIII_III.qs2"
)
qs2::qs_save(sce4.1,
  file = "results/sim_data/lm/sceIV_I.qs2"
)
qs2::qs_save(sce4.2,
  file = "results/sim_data/lm/sceIV_II.qs2"
)
qs2::qs_save(sce4.3,
  file = "results/sim_data/lm/sceIV_III.qs2"
)
