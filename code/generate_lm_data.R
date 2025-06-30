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
    X0k %*% beta0[, k] + matrix(rnorm(n0, sd = sqrt(sg0[k])), nrow = n0, ncol = 1)
  })
  # Generate current data
  n <- nrow(X)
  y <- X %*% beta + matrix(rnorm(n, sd = sqrt(sg)), nrow = n, ncol = 1)
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
# Generate perfect scenario 1
#-------------------------------------------------------------------------------

# define historical parameters
K <- 3
p0 <- 3
n0 <- 50

# define current parameters
n <- n0
p <- p0

# set number of datasets to generate
N <- 100

# set inputs to generate data
beta0 <- matrix(rep(c(1, -0.5, 0.5), K), nrow = p0, ncol = K)
sg0 <- rep(0.5, K)
inputs0 <- list(
  beta0 = beta0,
  sg0 = sg0,
  X0 = lapply(1:K, function(k) {
    cbind(1, matrix(rnorm(n0 * (p0-1)), nrow = n0, ncol = p0-1))
  }),
  beta = beta0[, 1],  # Assuming beta is the same as beta0 for simplicity
  sg = sg0[1],  # Assuming sg is the same as sg0 for simplicity
  X = cbind(1, matrix(rnorm(n * (p-1)), nrow = n, ncol = p-1))
)


inputs <- rep(list(inputs0), N)

# define cluster for parallel processing
cl <- makeCluster(detectCores() - 1)  # or 4, etc.
clusterExport(cl, varlist = c("flatten_data", "generate_lm"))  # Export function

# Generate datasets in parallel
sce1_data <- parLapply(cl, inputs, function(params) {
  generate_lm(params$X0, params$beta0, params$sg0,
              params$X, params$beta, params$sg)
})

# Stop the cluster
stopCluster(cl)


# Save the generated data
save(sce1_data, file = "results/sim_data/sce1_lm_data.RData")

#-------------------------------------------------------------------------------
# Generate perfect scenario 2
#-------------------------------------------------------------------------------

# define historical parameters
K <- 3
p0 <- 3
n0 <- 50

# define current parameters
n <- n0
p <- p0

# set number of datasets to generate
N <- 100

# set inputs to generate data
beta0 <- matrix(rep(c(1, -0.5, 0.5), K), nrow = p0, ncol = K)
sg0 <- rep(1, K)
inputs0 <- list(
  beta0 = beta0,
  sg0 = sg0,
  X0 = lapply(1:K, function(k) {
    cbind(1, matrix(rnorm(n0 * (p0-1)), nrow = n0, ncol = p0-1))
  }),
  beta = beta0[, 1],  # Assuming beta is the same as beta0 for simplicity
  sg = sg0[1],  # Assuming sg is the same as sg0 for simplicity
  X = cbind(1, matrix(rnorm(n * (p-1)), nrow = n, ncol = p-1))
)


inputs <- rep(list(inputs0), N)

# define cluster for parallel processing
cl <- makeCluster(detectCores() - 1)  # or 4, etc.
clusterExport(cl, varlist = c("flatten_data", "generate_lm"))  # Export function

# Generate datasets in parallel
sce2_data <- parLapply(cl, inputs, function(params) {
  generate_lm(params$X0, params$beta0, params$sg0,
              params$X, params$beta, params$sg)
})

# Stop the cluster
stopCluster(cl)


# Save the generated data
save(sce2_data, file = "results/sim_data/sce2_lm_data.RData")
