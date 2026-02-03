library(parallel)

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


# function to generate data
generate_lm <- function(X0, beta0, sg0, X, beta, sg) {
  # Generate historical data
  n0 <- nrow(X0[[1]])
  y0 <- lapply(1:length(X0), function(k) {
    X0k <- X0[[k]]
    X0k %*% beta0[, k] + matrix(rnorm(n0, sd = sg0[k]), nrow = n0, ncol = 1)
  })
  # Generate current data
  n <- nrow(X)
  y <- X %*% beta + matrix(rnorm(n, sd = sg), nrow = n, ncol = 1)
  # Flatten the data
  flattened_X0 <- flatten_data(X0)
  flattened_y0 <- flatten_data(y0)
  X0_flat <- flattened_X0$data_flat
  y0_flat <- flattened_y0$data_flat
  dims_X0 <- flattened_X0$dims_data
  dims_y0 <- flattened_y0$dims_data
  startid_X0 <- flattened_X0$starting_index
  startid_y0 <- flattened_y0$starting_index
  return(list(X0 = X0,
              X0_flat = X0_flat,
              dims_X0 = dims_X0,
              startid_X0 = startid_X0,
              y0 = y0,
              y0_flat = y0_flat,
              dims_y0 = dims_y0,
              startid_y0 = startid_y0,
              beta0 = beta0,
              sg0 = sg0,
              X = X, 
              y = y,
              beta = beta,
              sg = sg))
}

#-------------------------------------------------------------------------------
# Generate scenarios
#-------------------------------------------------------------------------------

# define historical parameters
K <- 3
p0 <- 3
n0 <- 100

# define current parameters
n <- n0
p <- p0

# set number of datasets to generate
N <- 200

betastar <- c(1, -0.5, 0.5)
sgstar <- 1
cs <- list(# Scenario I
           c(0.8, 0.9, 1),
           c(0.7, 0.8, 0.9),
           c(0.1, 0.3, 0.5),
           # Scenario II
           c(0.8, 1, 0.9),
           c(0.7, 0.9, 0.8),
           c(0.1, 0.5, 0.3),
           # Scenario III
           c(1, 0.9, 0.8),
           c(0.9, 0.8, 0.7),
           c(0.5, 0.3, 0.1),
           # Scenario IV
           c(1, 1, 1),
           c(0.6, 0.6, 0.6),
           c(0.2, 0.2, 0.2)
           )

generate_sce <- function(i) {
  c <- cs[[i]]
  # set inputs to generate data
  inputs0 <- list(
    beta0 = matrix(rep(betastar, K), nrow = p0, ncol = K) %*% diag(c),
    sg0 = as.vector(rep(sgstar, K) %*% diag(c)),
    X0 = lapply(1:K, function(k) {
      cbind(1, matrix(rnorm(n0 * (p0-1)), nrow = n0, ncol = p0-1))
    }),
    beta = betastar,  # Assuming beta is the same as beta0 for simplicity
    sg = sgstar,  # Assuming sg is the same as sg0 for simplicity
    X = cbind(1, matrix(rnorm(n * (p-1)), nrow = n, ncol = p-1))
  )
  
  
  inputs <- rep(list(inputs0), N)
  
  # Generate datasets
  sce_data <- mclapply(inputs, function(params) {
    generate_lm(params$X0, params$beta0, params$sg0,
                params$X, params$beta, params$sg)
  })
  sce <- as.roman(ceiling(i/3))
  sce <- ifelse(i %% 3 == 1, paste0(sce,"_I"), ifelse(i %% 3 == 2, paste0(sce,"_II"), paste0(sce,"_III")))
  saveRDS(sce_data, file = paste0("results/sim_data/lm/sce_", sce, ".rds"))
  return(sce_data)
}