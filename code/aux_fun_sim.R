# function to simulate data and run the models
sample_sce_bin <- function(par_list, gamma_model, delta_model, post = 1){
  # get parameters
  n0 <- par_list$n0
  n <- par_list$n
  a <- par_list$a
  b <- par_list$b
  al <- par_list$al
  bl <- par_list$bl
  theta0 <- par_list$theta0
  theta <- par_list$theta
  alpha <- par_list$alpha
  K <- length(n0)
  z0 <- par_list$z0
  z <- par_list$z
  if (is.null(z0)){
    z0 <- unlist(lapply(1:K, function(i) rbinom(1, n0[i], theta0[i])))
  }
  if (is.null(z)){
    z <- rbinom(1, n, theta)
  }
  # define data list
  data_onpp <- list(a = a, 
                    b = b, 
                    K = K, 
                    n = n, 
                    z = z, 
                    n0 = n0, 
                    z0 = z0, 
                    alpha = alpha,
                    post = post,
                    seq = 0
  )
  data_onppseq <- list(a = a, 
                       b = b, 
                       K = K, 
                       n = n, 
                       z = z, 
                       n0 = n0, 
                       z0 = z0, 
                       alpha = alpha,
                       post = post,
                       seq = 1
  )
  data_npp <- list(a = a, 
                   b = b, 
                   K = K, 
                   n = n, 
                   z = z, 
                   n0 = n0, 
                   z0 = z0, 
                   al = al,
                   bl = bl,
                   post = post,
                   seq = 0
  )
  data_nppseq <- list(a = a, 
                      b = b, 
                      K = K, 
                      n = n, 
                      z = z, 
                      n0 = n0, 
                      z0 = z0, 
                      al = al,
                      bl = bl,
                      post = post,
                      seq = 1
  )
  # sample from the models
  sample_delta_onpp <- gamma_model$sample(data = data_onpp, 
                                              chains = 4, 
                                              # parallel_chains = 4, 
                                              iter_warmup = 2000, 
                                              iter_sampling = 2000,
                                              adapt_delta = 0.98,
                                              refresh = 0
  )
  sample_delta_onppseq <- gamma_model$sample(data = data_onppseq, 
                                                 chains = 4, 
                                                 # parallel_chains = 4, 
                                                 iter_warmup = 2000, 
                                                 iter_sampling = 2000,
                                                 adapt_delta = 0.98,
                                                 refresh = 0
  )
  sample_delta_npp <- delta_model$sample(data = data_npp, 
                                             chains = 4, 
                                             # parallel_chains = 4, 
                                             iter_warmup = 2000, 
                                             iter_sampling = 2000,
                                             adapt_delta = 0.98,
                                             refresh = 0
  )
  sample_delta_nppseq <- delta_model$sample(data = data_nppseq, 
                                                chains = 4, 
                                                # parallel_chains = 4, 
                                                iter_warmup = 2000, 
                                                iter_sampling = 2000,
                                                adapt_delta = 0.98,
                                                refresh = 0
  )
  # # save output files
  # save_dir <- "results/cmdstan_outputs"
  # if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
  # sample_delta_npp$save_output_files(dir = save_dir)
  # sample_delta_nppseq$save_output_files(dir = save_dir)
  # sample_delta_onpp$save_output_files(dir = save_dir)
  # sample_delta_onppseq$save_output_files(dir = save_dir)
  
  # get diagnostics
  diagnostics_onpp <- sample_delta_onpp$sampler_diagnostics()
  diagnostics_onppseq <- sample_delta_onppseq$sampler_diagnostics()
  diagnostics_npp <- sample_delta_npp$sampler_diagnostics()
  diagnostics_nppseq <- sample_delta_nppseq$sampler_diagnostics()
  
  # get divergences
  divergences_onpp <- sum((diagnostics_onpp[, , "divergent__"] %>% as_draws_df())$divergent__)
  divergences_onppseq <- sum((diagnostics_onppseq[, , "divergent__"] %>% as_draws_df())$divergent__)
  divergences_npp <- sum((diagnostics_npp[, , "divergent__"] %>% as_draws_df())$divergent__)
  divergences_nppseq <- sum((diagnostics_nppseq[, , "divergent__"] %>% as_draws_df())$divergent__)
  
  # get draws
  draws_delta_onpp <- sample_delta_onpp$draws(variables = "delta") %>% as_draws_matrix()
  draws_delta_onppseq <- sample_delta_onppseq$draws(variables = "delta") %>% as_draws_matrix()
  draws_delta_npp <- sample_delta_npp$draws(variables = "delta") %>% as_draws_matrix()
  draws_delta_nppseq <- sample_delta_nppseq$draws(variables = "delta") %>% as_draws_matrix()
  
  theta_npp <- rbeta(nrow(draws_delta_npp),
                     z + a + draws_delta_npp %*% z0, 
                     n - z + b + draws_delta_npp %*% (n0 - z0)
  )
  theta_nppseq <- rbeta(nrow(draws_delta_nppseq),
                        z + K*(a-1) + 1 + draws_delta_nppseq %*% z0, 
                        n - z + K*(b-1) + 1 + draws_delta_nppseq %*% (n0 - z0)
  )
  theta_onpp <- rbeta(nrow(draws_delta_onpp),
                      z + a + draws_delta_onpp %*% z0, 
                      n - z + b + draws_delta_onpp %*% (n0 - z0)
  )
  theta_onppseq <- rbeta(nrow(draws_delta_onppseq),
                         z + K*(a-1) + 1 + draws_delta_onppseq %*% z0, 
                         n - z + K*(b-1) + 1 + draws_delta_onppseq %*% (n0 - z0)
  )
  
  return(list(hattheta_npp = mean(theta_npp),
              hattheta_nppseq = mean(theta_nppseq),
              hattheta_onpp = mean(theta_onpp),
              hattheta_onppseq = mean(theta_onppseq),
              hatdelta_npp = colMeans(draws_delta_npp),
              hatdelta_nppseq = colMeans(draws_delta_nppseq),
              hatdelta_onpp = colMeans(draws_delta_onpp),
              hatdelta_onppseq = colMeans(draws_delta_onppseq),
              divergences = list(divergences_npp = divergences_npp,
                                 divergences_nppseq = divergences_nppseq,
                                 divergences_onpp = divergences_onpp,
                                 divergences_onppseq = divergences_onppseq),
              delta_npp = draws_delta_npp,
              delta_nppseq = draws_delta_nppseq,
              delta_onpp = draws_delta_onpp,
              delta_onppseq = draws_delta_onppseq,
              theta_npp = theta_npp,
              theta_nppseq = theta_nppseq,
              theta_onpp = theta_onpp,
              theta_onppseq = theta_onppseq
  ))
}

# function to simulate data and run the models
sample_sce_poi <- function(par_list, gamma_model, delta_model){
  # get parameters
  t0 <- par_list$t0
  t <- par_list$t
  a <- par_list$a
  b <- par_list$b
  al <- par_list$al
  bl <- par_list$bl
  theta0 <- par_list$theta0
  theta <- par_list$theta
  alpha <- par_list$alpha
  K <- length(t0)
  z0 <- par_list$z0
  z <- par_list$z
  if (is.null(z0)){
    z0 <- unlist(lapply(1:K, function(i) rpois(1, t0[i]*theta0[i])))
  }
  if (is.null(z)){
    z <- rpois(1, t*theta)
  }
  # define data list
  data_onpp <- list(a = a, 
                    b = b, 
                    K = K, 
                    t = t, 
                    z = z, 
                    t0 = t0, 
                    z0 = z0, 
                    alpha = alpha,
                    post = 1,
                    seq = 0
  )
  data_onppseq <- list(a = a, 
                       b = b, 
                       K = K, 
                       t = t, 
                       z = z, 
                       t0 = t0, 
                       z0 = z0, 
                       alpha = alpha,
                       post = 1,
                       seq = 1
  )
  data_npp <- list(a = a, 
                   b = b, 
                   K = K, 
                   t = t, 
                   z = z, 
                   t0 = t0, 
                   z0 = z0, 
                   al = al,
                   bl = bl,
                   post = 1,
                   seq = 0
  )
  data_nppseq <- list(a = a, 
                      b = b, 
                      K = K, 
                      t = t, 
                      z = z, 
                      t0 = t0, 
                      z0 = z0, 
                      al = al,
                      bl = bl,
                      post = 1,
                      seq = 1
  )
  # sample from the models
  sample_delta_onpp <- gamma_model$sample(data = data_onpp, 
                                          chains = 4, 
                                          iter_warmup = 2000, 
                                          iter_sampling = 2000,
                                          adapt_delta = 0.98,
                                          refresh = 0
  )
  sample_delta_onppseq <- gamma_model$sample(data = data_onppseq, 
                                             chains = 4, 
                                             iter_warmup = 2000, 
                                             iter_sampling = 2000,
                                             adapt_delta = 0.98,
                                             refresh = 0
  )
  sample_delta_npp <- delta_model$sample(data = data_npp, 
                                         chains = 4, 
                                         iter_warmup = 2000, 
                                         iter_sampling = 2000,
                                         adapt_delta = 0.98,
                                         refresh = 0
  )
  sample_delta_nppseq <- delta_model$sample(data = data_nppseq, 
                                            chains = 4,  
                                            iter_warmup = 2000, 
                                            iter_sampling = 2000,
                                            adapt_delta = 0.98,
                                            refresh = 0
  )
  # # save output files
  # save_dir <- "results/cmdstan_outputs"
  # if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
  # sample_delta_npp$save_output_files(dir = save_dir)
  # sample_delta_nppseq$save_output_files(dir = save_dir)
  # sample_delta_onpp$save_output_files(dir = save_dir)
  # sample_delta_onppseq$save_output_files(dir = save_dir)
  
  # get diagnostics
  diagnostics_onpp <- sample_delta_onpp$sampler_diagnostics()
  diagnostics_onppseq <- sample_delta_onppseq$sampler_diagnostics()
  diagnostics_npp <- sample_delta_npp$sampler_diagnostics()
  diagnostics_nppseq <- sample_delta_nppseq$sampler_diagnostics()
  
  # get divergences
  divergences_onpp <- sum((diagnostics_onpp[, , "divergent__"] %>% as_draws_df())$divergent__)
  divergences_onppseq <- sum((diagnostics_onppseq[, , "divergent__"] %>% as_draws_df())$divergent__)
  divergences_npp <- sum((diagnostics_npp[, , "divergent__"] %>% as_draws_df())$divergent__)
  divergences_nppseq <- sum((diagnostics_nppseq[, , "divergent__"] %>% as_draws_df())$divergent__)
  
  # get draws
  draws_delta_onpp <- sample_delta_onpp$draws(variables = "delta") %>% as_draws_matrix()
  draws_delta_onppseq <- sample_delta_onppseq$draws(variables = "delta") %>% as_draws_matrix()
  draws_delta_npp <- sample_delta_npp$draws(variables = "delta") %>% as_draws_matrix()
  draws_delta_nppseq <- sample_delta_nppseq$draws(variables = "delta") %>% as_draws_matrix()
  
  theta_npp <- rgamma(nrow(draws_delta_npp),
                      z + a + draws_delta_npp %*% z0, 
                      b + draws_delta_npp %*% t0 + t
  )
  theta_nppseq <- rgamma(nrow(draws_delta_nppseq),
                         z + K*(a-1) + 1 + draws_delta_nppseq %*% z0, 
                         K*b + draws_delta_nppseq %*% t0 + t
  )
  theta_onpp <- rgamma(nrow(draws_delta_onpp),
                       z + a + draws_delta_onpp %*% z0, 
                       b + draws_delta_onpp %*% t0 + t
  )
  theta_onppseq <- rgamma(nrow(draws_delta_onppseq),
                          z + K*(a-1) + 1 + draws_delta_onppseq %*% z0, 
                          K*b + draws_delta_onppseq %*% t0 + t
  )
  
  return(list(hattheta_npp = mean(theta_npp),
              hattheta_nppseq = mean(theta_nppseq),
              hattheta_onpp = mean(theta_onpp),
              hattheta_onppseq = mean(theta_onppseq),
              hatdelta_npp = colMeans(draws_delta_npp),
              hatdelta_nppseq = colMeans(draws_delta_nppseq),
              hatdelta_onpp = colMeans(draws_delta_onpp),
              hatdelta_onppseq = colMeans(draws_delta_onppseq),
              divergences = list(divergences_npp = divergences_npp,
                                 divergences_nppseq = divergences_nppseq,
                                 divergences_onpp = divergences_onpp,
                                 divergences_onppseq = divergences_onppseq),
              delta_npp = draws_delta_npp,
              delta_nppseq = draws_delta_nppseq,
              delta_onpp = draws_delta_onpp,
              delta_onppseq = draws_delta_onppseq,
              theta_npp = theta_npp,
              theta_nppseq = theta_nppseq,
              theta_onpp = theta_onpp,
              theta_onppseq = theta_onppseq
  ))
}

sample_sce_normal_fixed_var <- function(par_list, gamma_model, delta_model, post = 1){
  # get parameters
  n0 <- par_list$n0
  n <- par_list$n
  mu0 <- par_list$mu0
  sigma0 <- par_list$sigma0
  al <- par_list$al
  bl <- par_list$bl
  alpha <- par_list$alpha
  theta0 <- par_list$theta0
  theta <- par_list$theta
  sigmah <- par_list$sigmah
  sigma <- par_list$sigma
  
  K <- length(n0)
  y0 <- par_list$y0
  y <- par_list$y
  if (is.null(y0)){
    y0 <- do.call(c, lapply(1:K, function(i) 
                                  rnorm(n0[i], mean = theta0[i], sd = sigmah[i])))
    start_idx <- c(1, cumsum(n0) + 1)[1:K]
  }
  if (is.null(y)){
    y <- rnorm(n, mean = theta, sd = sigma)
  }
  # define data list
  data_onpp <- list(
    K = K,
    n0 = n0,
    n = n,
    mu0 = mu0,
    sigma0 = sigma0,
    sigmah = sigmah,
    sigma = sigma,
    y0 = y0,
    y = y,
    start_idx = start_idx,
    alpha = alpha,
    post = post,
    seq = 0
  )
  data_onppseq <- list(
    K = K,
    n0 = n0,
    n = n,
    mu0 = mu0,
    sigma0 = sigma0,
    sigmah = sigmah,
    sigma = sigma,
    y0 = y0,
    y = y,
    start_idx = start_idx,
    alpha = alpha,
    post = post,
    seq = 1
  )
  data_npp <- list(
    K = K,
    n0 = n0,
    n = n,
    mu0 = mu0,
    sigma0 = sigma0,
    sigmah = sigmah,
    sigma = sigma,
    y0 = y0,
    y = y,
    start_idx = start_idx,
    al = al,
    bl = bl,
    post = post,
    seq = 0
  )
  data_nppseq <- list(
    K = K,
    n0 = n0,
    n = n,
    mu0 = mu0,
    sigma0 = sigma0,
    sigmah = sigmah,
    sigma = sigma,
    y0 = y0,
    y = y,
    start_idx = start_idx,
    al = al,
    bl = bl,
    post = post,
    seq = 1
  )
  # sample from the models
  sample_delta_onpp <- gamma_model$sample(data = data_onpp, 
                                              chains = 4, 
                                              # parallel_chains = 4, 
                                              iter_warmup = 2000, 
                                              iter_sampling = 2000,
                                              adapt_delta = 0.9999999,
                                              refresh = 0
  )
  sample_delta_onppseq <- gamma_model$sample(data = data_onppseq, 
                                                 chains = 4, 
                                                 # parallel_chains = 4, 
                                                 iter_warmup = 2000, 
                                                 iter_sampling = 2000,
                                                 adapt_delta = 0.999999,
                                                 refresh = 0
  )
  sample_delta_npp <- delta_model$sample(data = data_npp, 
                                             chains = 4, 
                                             # parallel_chains = 4, 
                                             iter_warmup = 2000, 
                                             iter_sampling = 2000,
                                             adapt_delta = 0.999999,
                                             refresh = 0
  )
  sample_delta_nppseq <- delta_model$sample(data = data_nppseq, 
                                                chains = 4, 
                                                # parallel_chains = 4, 
                                                iter_warmup = 2000, 
                                                iter_sampling = 2000,
                                                adapt_delta = 0.999999,
                                                refresh = 0
  )
  # # save output files
  # save_dir <- "results/cmdstan_outputs"
  # if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
  # sample_delta_npp$save_output_files(dir = save_dir)
  # sample_delta_nppseq$save_output_files(dir = save_dir)
  # sample_delta_onpp$save_output_files(dir = save_dir)
  # sample_delta_onppseq$save_output_files(dir = save_dir)
  
  # get diagnostics
  diagnostics_onpp <- sample_delta_onpp$sampler_diagnostics()
  diagnostics_onppseq <- sample_delta_onppseq$sampler_diagnostics()
  diagnostics_npp <- sample_delta_npp$sampler_diagnostics()
  diagnostics_nppseq <- sample_delta_nppseq$sampler_diagnostics()
  
  # get divergences
  divergences_onpp <- sum((diagnostics_onpp[, , "divergent__"] %>% as_draws_df())$divergent__)
  divergences_onppseq <- sum((diagnostics_onppseq[, , "divergent__"] %>% as_draws_df())$divergent__)
  divergences_npp <- sum((diagnostics_npp[, , "divergent__"] %>% as_draws_df())$divergent__)
  divergences_nppseq <- sum((diagnostics_nppseq[, , "divergent__"] %>% as_draws_df())$divergent__)
  
  # get draws
  draws_delta_onpp <- sample_delta_onpp$draws(variables = "delta") %>% as_draws_matrix()
  draws_delta_onppseq <- sample_delta_onppseq$draws(variables = "delta") %>% as_draws_matrix()
  draws_delta_npp <- sample_delta_npp$draws(variables = "delta") %>% as_draws_matrix()
  draws_delta_nppseq <- sample_delta_nppseq$draws(variables = "delta") %>% as_draws_matrix()
  
  theta_npp <- sample_delta_npp$draws(variables = "theta") %>% as_draws_df() %>% pull(theta)
  theta_nppseq <- sample_delta_nppseq$draws(variables = "theta") %>% as_draws_df() %>% pull(theta)
  theta_onpp <- sample_delta_onpp$draws(variables = "theta") %>% as_draws_df() %>% pull(theta)
  theta_onppseq <- sample_delta_onppseq$draws(variables = "theta") %>% as_draws_df() %>% pull(theta)
  
  return(list(hattheta_npp = mean(theta_npp),
              hattheta_nppseq = mean(theta_nppseq),
              hattheta_onpp = mean(theta_onpp),
              hattheta_onppseq = mean(theta_onppseq),
              hatdelta_npp = colMeans(draws_delta_npp),
              hatdelta_nppseq = colMeans(draws_delta_nppseq),
              hatdelta_onpp = colMeans(draws_delta_onpp),
              hatdelta_onppseq = colMeans(draws_delta_onppseq),
              divergences = list(divergences_npp = divergences_npp,
                                 divergences_nppseq = divergences_nppseq,
                                 divergences_onpp = divergences_onpp,
                                 divergences_onppseq = divergences_onppseq),
              delta_npp = draws_delta_npp,
              delta_nppseq = draws_delta_nppseq,
              delta_onpp = draws_delta_onpp,
              delta_onppseq = draws_delta_onppseq,
              theta_npp = theta_npp,
              theta_nppseq = theta_nppseq,
              theta_onpp = theta_onpp,
              theta_onppseq = theta_onppseq
  ))
}

sample_sce_lm <- function(par_list, data, gamma_model, delta_model, post = 1){
  X0 <- data$X0_flat
  startid_X0 <- data$startid_X0
  n0 <- as.vector(data$dims_y0[,1]) 
  y0 <- data$y0_flat
  startid_y0 <- data$startid_y0
  K <- length(startid_X0)
  X <- data$X
  p <- ncol(X)
  n <- nrow(X)
  len_X0 <- length(X0)
  len_y0 <- length(y0)
  y <- as.vector(data$y)
  # get parameters
  mu0 <- par_list$mu0
  V0 <- par_list$V0
  a <- par_list$a
  b <- par_list$b
  tilde_a <- par_list$tilde_a
  tilde_b <- par_list$tilde_b
  alpha <- par_list$alpha

  # define data list
  data_onpp <- list(
    K = K,
    n0 = n0,
    n = n,
    p = p,
    len_X0 = len_X0,
    len_y0 = len_y0,
    X0 = X0,
    startid_X0 = startid_X0,
    y0 = y0,
    startid_y0 = startid_y0,
    X = X,
    y = y,
    a = a,
    b = b,
    V0 = V0,
    mu0 = mu0,
    alpha = alpha,
    post = post,
    seq = 0
  )
  data_onppseq <- list(
    K = K,
    n0 = n0,
    n = n,
    p = p,
    len_X0 = len_X0,
    len_y0 = len_y0,
    X0 = X0,
    startid_X0 = startid_X0,
    y0 = y0,
    startid_y0 = startid_y0,
    X = X,
    y = y,
    a = a,
    b = b,
    V0 = V0,
    mu0 = mu0,
    alpha = alpha,
    post = post,
    seq = 1
  )
  data_npp <- list(
    K = K,
    n0 = n0,
    n = n,
    p = p,
    len_X0 = len_X0,
    len_y0 = len_y0,
    X0 = X0,
    startid_X0 = startid_X0,
    y0 = y0,
    startid_y0 = startid_y0,
    X = X,
    y = y,
    a = a,
    b = b,
    V0 = V0,
    mu0 = mu0,
    tilde_a = tilde_a,
    tilde_b = tilde_b,
    post = post,
    seq = 0
  )
  data_nppseq <- list(
    K = K,
    n0 = n0,
    n = n,
    p = p,
    len_X0 = len_X0,
    len_y0 = len_y0,
    X0 = X0,
    startid_X0 = startid_X0,
    y0 = y0,
    startid_y0 = startid_y0,
    X = X,
    y = y,
    a = a,
    b = b,
    V0 = V0,
    mu0 = mu0,
    tilde_a = tilde_a,
    tilde_b = tilde_b,
    post = post,
    seq = 1
  )
  # sample from the models
  sample_delta_onpp <- gamma_model$sample(data = data_onpp, 
                                              chains = 4, 
                                              # parallel_chains = 4, 
                                              iter_warmup = 2000, 
                                              iter_sampling = 2000,
                                              # adapt_delta = 0.9999999,
                                              refresh = 0
  )
  sample_delta_onppseq <- gamma_model$sample(data = data_onppseq, 
                                                 chains = 4, 
                                                 # parallel_chains = 4, 
                                                 iter_warmup = 2000, 
                                                 iter_sampling = 2000,
                                                #  adapt_delta = 0.999999,
                                                 refresh = 0
  )
  sample_delta_npp <- delta_model$sample(data = data_npp, 
                                             chains = 4, 
                                             # parallel_chains = 4, 
                                             iter_warmup = 2000, 
                                             iter_sampling = 2000,
                                            #  adapt_delta = 0.999999,
                                             refresh = 0
  )
  sample_delta_nppseq <- delta_model$sample(data = data_nppseq, 
                                                chains = 4, 
                                                # parallel_chains = 4, 
                                                iter_warmup = 2000, 
                                                iter_sampling = 2000,
                                                # adapt_delta = 0.999999,
                                                refresh = 0
  )
  # # save output files
  # save_dir <- "results/cmdstan_outputs"
  # if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
  # sample_delta_npp$save_output_files(dir = save_dir)
  # sample_delta_nppseq$save_output_files(dir = save_dir)
  # sample_delta_onpp$save_output_files(dir = save_dir)
  # sample_delta_onppseq$save_output_files(dir = save_dir)
  
  # get diagnostics
  diagnostics_onpp <- sample_delta_onpp$sampler_diagnostics()
  diagnostics_onppseq <- sample_delta_onppseq$sampler_diagnostics()
  diagnostics_npp <- sample_delta_npp$sampler_diagnostics()
  diagnostics_nppseq <- sample_delta_nppseq$sampler_diagnostics()
  
  # get divergences
  divergences_onpp <- sum((diagnostics_onpp[, , "divergent__"] %>% as_draws_df())$divergent__)
  divergences_onppseq <- sum((diagnostics_onppseq[, , "divergent__"] %>% as_draws_df())$divergent__)
  divergences_npp <- sum((diagnostics_npp[, , "divergent__"] %>% as_draws_df())$divergent__)
  divergences_nppseq <- sum((diagnostics_nppseq[, , "divergent__"] %>% as_draws_df())$divergent__)
  
  # get draws
  draws_delta_onpp <- sample_delta_onpp$draws(variables = "delta") %>% as_draws_matrix()
  draws_delta_onppseq <- sample_delta_onppseq$draws(variables = "delta") %>% as_draws_matrix()
  draws_delta_npp <- sample_delta_npp$draws(variables = "delta") %>% as_draws_matrix()
  draws_delta_nppseq <- sample_delta_nppseq$draws(variables = "delta") %>% as_draws_matrix()
  
  theta_npp <- sample_delta_npp$draws(variables = "theta") %>% as_draws_matrix() 
  theta_nppseq <- sample_delta_nppseq$draws(variables = "theta") %>% as_draws_matrix()
  theta_onpp <- sample_delta_onpp$draws(variables = "theta") %>% as_draws_matrix()
  theta_onppseq <- sample_delta_onppseq$draws(variables = "theta") %>% as_draws_matrix()
  
  return(list(hattheta_npp = colMeans(theta_npp),
              hattheta_nppseq = colMeans(theta_nppseq),
              hattheta_onpp = colMeans(theta_onpp),
              hattheta_onppseq = colMeans(theta_onppseq),
              hatdelta_npp = colMeans(draws_delta_npp),
              hatdelta_nppseq = colMeans(draws_delta_nppseq),
              hatdelta_onpp = colMeans(draws_delta_onpp),
              hatdelta_onppseq = colMeans(draws_delta_onppseq),
              divergences = list(divergences_npp = divergences_npp,
                                 divergences_nppseq = divergences_nppseq,
                                 divergences_onpp = divergences_onpp,
                                 divergences_onppseq = divergences_onppseq),
              delta_npp = draws_delta_npp,
              delta_nppseq = draws_delta_nppseq,
              delta_onpp = draws_delta_onpp,
              delta_onppseq = draws_delta_onppseq,
              theta_npp = theta_npp,
              theta_nppseq = theta_nppseq,
              theta_onpp = theta_onpp,
              theta_onppseq = theta_onppseq
  ))
}

# function to run simulations in parallel
sim_sce <- function(model, num_cores, num_sim, sce, gamma_model, delta_model){
  # Start a cluster with the desired number of cores
  cl <- makeCluster(num_cores)
  # Load the required packages
  clusterEvalQ(cl, {
    library(cmdstanr)
    library(MCMCpack)
    library(dplyr)
    library(tidyverse)
    library(posterior)
  })
  # Export the models to the cluster
  clusterExport(cl, 
                varlist = c("gamma_model", 
                            "delta_model"),
                envir = environment()
  )
  
  # Run simulations in parallel
  if (model == "bin"){
    sample_sce <- sample_sce_bin
    rep_sce <- replicate(num_sim, sce, simplify = FALSE)
  } else if (model == "poi"){
    sample_sce <- sample_sce_poi
    rep_sce <- replicate(num_sim, sce, simplify = FALSE)
  } else if (model == "normal_fixed_var"){
    sample_sce <- sample_sce_normal_fixed_var
    rep_sce <- replicate(num_sim, sce, simplify = FALSE)
  } else if (model == "lm"){
    sample_sce <- sample_sce_lm
    rep_sce <- sce$data
  } else {
    stop("Invalid model type. Choose from 'bin', 'poi', 'normal_fixed_var', or 'lm'.")
  }

  pbapply::pboptions(type = "timer")
  if (model == "lm") {
    results <- pbapply::pblapply( 
                       rep_sce, 
                       par_list = sce$par,
                       sample_sce,
                       gamma_model = gamma_model,
                       delta_model = delta_model,
                       cl = cl
                      )
  } else {
    results <- pbapply::pblapply( 
                       rep_sce, 
                       sample_sce,
                       gamma_model = gamma_model,
                       delta_model = delta_model,
                       cl = cl
                      )
  }
  
  
  # Stop the cluster after computation
  stopCluster(cl)
  
  # combine results
  hattheta <- data.frame(
    hattheta_npp = sapply(results, function(x) x$hattheta_npp),
    hattheta_nppseq = sapply(results, function(x) x$hattheta_nppseq),
    hattheta_onpp = sapply(results, function(x) x$hattheta_onpp),
    hattheta_onppseq = sapply(results, function(x) x$hattheta_onppseq)
  )
  
  divergences <- data.frame(divergences_npp = sum(sapply(results, function(x) x$divergences$divergences_npp)),
                            divergences_nppseq = sum(sapply(results, function(x) x$divergences$divergences_nppseq)),
                            divergences_onpp = sum(sapply(results, function(x) x$divergences$divergences_onpp)),
                            divergences_onppseq = sum(sapply(results, function(x) x$divergences$divergences_onppseq))
  )
  
  hatdelta <- list()
  for (i in 1:length(results[[1]]$hatdelta_npp)){
    hatdelta[[i]] <- data.frame(
      hatdelta_npp = sapply(results, function(x) x$hatdelta_npp[i]),
      hatdelta_nppseq = sapply(results, function(x) x$hatdelta_nppseq[i]),
      hatdelta_onpp = sapply(results, function(x) x$hatdelta_onpp[i]),
      hatdelta_onppseq = sapply(results, function(x) x$hatdelta_onppseq[i])
    )
  }
  
  return(list(hattheta = hattheta, 
              hatdelta = hatdelta,
              divergences = divergences,
              delta = lapply(results, function(x) list(delta_npp = x$delta_npp,
                                                        delta_nppseq = x$delta_nppseq,
                                                        delta_onpp = x$delta_onpp,
                                                        delta_onppseq = x$delta_onppseq)
              ),
              theta = lapply(results, function(x) list(theta_npp = x$theta_npp,
                                                        theta_nppseq = x$theta_nppseq,
                                                        theta_onpp = x$theta_onpp,
                                                        theta_onppseq = x$theta_onppseq)
              )
  ))
}

# Extract draws from the samples
get_hat_par <- function(samples) {
  hat_delta_npp <- do.call(rbind, 
                           lapply(samples, 
                                  function(x) colMeans(x$draws_delta_npp)))
  hat_delta_nppseq <- do.call(rbind,
                              lapply(samples, 
                                     function(x) colMeans(x$draws_delta_nppseq)))
  hat_delta_onpp <- do.call(rbind,
                            lapply(samples, 
                                   function(x) colMeans(x$draws_delta_onpp)))
  hat_delta_onppseq <- do.call(rbind,
                               lapply(samples, 
                                      function(x) colMeans(x$draws_delta_onppseq)))
  hat_beta_npp <- do.call(rbind,
                          lapply(samples, 
                                 function(x) colMeans(x$draws_beta_npp)))
  hat_beta_nppseq <- do.call(rbind,
                             lapply(samples, 
                                    function(x) colMeans(x$draws_beta_nppseq)))
  hat_beta_onpp <- do.call(rbind,
                           lapply(samples, 
                                  function(x) colMeans(x$draws_beta_onpp)))
  hat_beta_onppseq <- do.call(rbind,
                              lapply(samples, 
                                     function(x) colMeans(x$draws_beta_onppseq)))
  hat_sigma_npp <- do.call(rbind,
                           lapply(samples, 
                                  function(x) colMeans(x$draws_sigma_npp)))
  hat_sigma_nppseq <- do.call(rbind,
                              lapply(samples, 
                                     function(x) colMeans(x$draws_sigma_nppseq)))
  hat_sigma_onpp <- do.call(rbind,
                            lapply(samples, 
                                   function(x) colMeans(x$draws_sigma_onpp)))
  hat_sigma_onppseq <- do.call(rbind,
                               lapply(samples, 
                                      function(x) colMeans(x$draws_sigma_onppseq)))
  hat_delta <- lapply(seq_len(ncol(hat_delta_npp)), function(i) {
    cbind(hat_delta_npp[,i],
          hat_delta_nppseq[,i],
          hat_delta_onpp[,i],
          hat_delta_onppseq[,i])
  })
  hat_theta <- lapply(1:(ncol(hat_beta_npp) + 1), function(i) {
    if (i <= ncol(hat_beta_npp)) {
      cbind(hat_beta_npp[,i],
            hat_beta_nppseq[,i],
            hat_beta_onpp[,i],
            hat_beta_onppseq[,i])
    } else {
      cbind(hat_sigma_npp,
            hat_sigma_nppseq,
            hat_sigma_onpp,
            hat_sigma_onppseq)
    }
    
  })
  return(list(hatdelta = hat_delta,
              hattheta = hat_theta))
}

# function to plot boxplots
plot_boxplot <- function(draws, variable) {
  # create data frame
  nc1 <- length(draws[,1])
  nc2 <- length(draws[,2])
  nc3 <- length(draws[,3])
  nc4 <- length(draws[,4])
  
  c <- data.frame(x=c(rep("NPP1",nc1),
                      rep("NPP2",nc2), 
                      rep("ONPP1",nc3),
                      rep("ONPP2",nc4)), 
                  y=c(draws[,1],
                      draws[,2],
                      draws[,3],
                      draws[,4])
  )
  
  c$x <- factor(c$x,levels=c('NPP1','NPP2','ONPP1','ONPP2'),ordered=TRUE)
  
  # plot
  box <- ggplot(c, aes(x=x, y=y, fill=x)) +
    theme_bw()+labs(x = NULL) + 
    labs(y = variable)+
    stat_boxplot(geom ='errorbar', coef = 0.5, width = 0.25)+
    geom_boxplot(outlier.colour="black", 
                 width = 0.4, outlier.size = 0.3, coef = 0.5)+
    stat_summary(fun = mean, colour="yellow3", geom="point", 
                 shape=18, size=3, show.legend = FALSE)+
    theme(legend.position="none",
          axis.title.y=element_text(angle= -270,  face='bold', size=11))+
    scale_fill_manual(values=c("#999999","red2", "green4", "dodgerblue1"))+
    scale_x_discrete(labels=c(expression(NPP), expression(NPP[seq]), 
                              expression(ONPP), expression(ONPP[seq]))
    )
}

# function to compute BCI
compute_bci <- function(list_draws, alpha){
  bci_npp <- lapply(list_draws, function(x) quantile(x$theta_npp, 
                                                     probs = c((1-alpha)/2, (1+alpha)/2)))
  bci_nppseq <- lapply(list_draws, function(x) quantile(x$theta_nppseq, 
                                                        probs = c((1-alpha)/2, (1+alpha)/2)))
  bci_onpp <- lapply(list_draws, function(x) quantile(x$theta_onpp,
                                                      probs = c((1-alpha)/2, (1+alpha)/2)))
  bci_onppseq <- lapply(list_draws, function(x) quantile(x$theta_onppseq,
                                                         probs = c((1-alpha)/2, (1+alpha)/2)))
  return(list(bci_npp = bci_npp,
              bci_nppseq = bci_nppseq,
              bci_onpp = bci_onpp,
              bci_onppseq = bci_onppseq
  ))
}

plot_sce_bin <- function(j) {
  sim <- sim_sces[[j]]
  sce <- as.roman(ceiling(j/3))
  sce <- ifelse(j %% 3 == 1, paste0(sce,".I"), ifelse(j %% 3 == 2, paste0(sce,".II"), paste0(sce,".III")))
  plots_delta <- lapply(seq_along(sim$hatdelta), function(i) {
    plot <- plot_boxplot(sim$hatdelta[[i]], bquote(hat(delta)[.(i)])) + 
      coord_flip() +
      ggtitle(paste0("Scenario ", sce))
      
    if (i>1) {
      plot <- plot + scale_x_discrete(labels = NULL) + ggtitle(NULL)
    }
    return(plot)
  })
  plot_theta <- plot_boxplot(sim$hattheta, expression(hat(theta))) + 
    geom_hline(yintercept = true_value, linetype = "dotted", color = "red", size = 1) +
    coord_flip()
  # Combine plots
  plot <- Reduce(`+`, plots_delta) + plot_theta + plot_layout(ncol = length(plots_delta) + 1) + 
    scale_x_discrete(labels = NULL)
  return(plot)
}

# function to compute coverage
compute_coverage <- function(bcis, true_theta) {
  coverage_npp <- sum(sapply(bcis$bci_npp, function(x) x[1] <= true_theta & true_theta <= x[2]))/
    length(bcis$bci_npp)
  coverage_nppseq <- sum(sapply(bcis$bci_nppseq, function(x) x[1] <= true_theta & true_theta <= x[2]))/
    length(bcis$bci_nppseq)
  coverage_onpp <- sum(sapply(bcis$bci_onpp, function(x) x[1] <= true_theta & true_theta <= x[2]))/
    length(bcis$bci_onpp)
  coverage_onppseq <- sum(sapply(bcis$bci_onppseq, function(x) x[1] <= true_theta & true_theta <= x[2]))/
    length(bcis$bci_onppseq)
  return(list(coverage_npp = coverage_npp,
              coverage_nppseq = coverage_nppseq,
              coverage_onpp = coverage_onpp,
              coverage_onppseq = coverage_onppseq
  ))
}

compute_avg_len <- function(bcis) {
  bcis_npp <- do.call(rbind,bcis$bci_npp)
  bcis_nppseq <- do.call(rbind,bcis$bci_nppseq)
  bcis_onpp <- do.call(rbind,bcis$bci_onpp)
  bcis_onppseq <- do.call(rbind,bcis$bci_onppseq)
  return(list(
    npp = mean(bcis_npp[,2] - bcis_npp[,1]),
    nppseq = mean(bcis_nppseq[,2] - bcis_nppseq[,1]),
    onpp = mean(bcis_onpp[,2] - bcis_onpp[,1]),
    onppseq = mean(bcis_onppseq[,2] - bcis_onppseq[,1])
  ))
}

compute_wis <- function(list_draws, quantile_level, true_value) {
  N <- length(list_draws)
  q_npp <- do.call(rbind,
                   lapply(list_draws, function(x) quantile(x$theta_npp, 
                                                           probs = quantile_level))
  )
  q_nppseq <- do.call(rbind,
                      lapply(list_draws, function(x) quantile(x$theta_nppseq, 
                                                              probs = quantile_level))
  )
  q_onpp <- do.call(rbind,
                    lapply(list_draws, function(x) quantile(x$theta_onpp,
                                                            probs = quantile_level))
  )
  q_onppseq <- do.call(rbind,
                       lapply(list_draws, function(x) quantile(x$theta_onppseq,
                                                               probs = quantile_level))
  )
  wis_npp <- wis(rep(true_value, N), q_npp, quantile_level)
  wis_nppseq <- wis(rep(true_value, N), q_nppseq, quantile_level)
  wis_onpp <- wis(rep(true_value, N), q_onpp, quantile_level)
  wis_onppseq <- wis(rep(true_value, N), q_onppseq, quantile_level)
  return(list(hat_wis_npp = mean(wis_npp),
              hat_wis_nppseq = mean(wis_nppseq),
              hat_wis_onpp = mean(wis_onpp),
              hat_wis_onppseq = mean(wis_onppseq)
  ))
}

compute_wis_lm <- function(list_draws, quantile_level, true_value) {
  N <- length(list_draws)
  q_npp <- lapply(1:(ncol(list_draws[[1]]$draws_beta_npp) + 1), function(j) {
    if (j <= ncol(list_draws[[1]]$draws_beta_npp)) {
      do.call(rbind,
              lapply(list_draws, function(x) {
                quantile(x$draws_beta_npp[, j], probs = quantile_level)
              })
      )
    } else {
      do.call(rbind,
              lapply(list_draws, function(x) {
                quantile(x$draws_sigma_npp, probs = quantile_level)
              })
      )
    }
  })
  q_nppseq <- lapply(1:(ncol(list_draws[[1]]$draws_beta_nppseq) + 1), function(j) {
    if (j <= ncol(list_draws[[1]]$draws_beta_nppseq)) {
      do.call(rbind,
              lapply(list_draws, function(x) {
                quantile(x$draws_beta_nppseq[, j], probs = quantile_level)
              })
      )
    } else {
      do.call(rbind,
              lapply(list_draws, function(x) {
                quantile(x$draws_sigma_nppseq, probs = quantile_level)
              })
      )
    }
  }
  )
  q_onpp <- lapply(1:(ncol(list_draws[[1]]$draws_beta_onpp) + 1), function(j) {
    if (j <= ncol(list_draws[[1]]$draws_beta_onpp)) {
      do.call(rbind,
              lapply(list_draws, function(x) {
                quantile(x$draws_beta_onpp[, j], probs = quantile_level)
              })
      )
    } else {
      do.call(rbind,
              lapply(list_draws, function(x) {
                quantile(x$draws_sigma_onpp, probs = quantile_level)
              })
      )
    }
  })
  q_onppseq <- lapply(1:(ncol(list_draws[[1]]$draws_beta_onppseq) + 1), function(j) {
    if (j <= ncol(list_draws[[1]]$draws_beta_onppseq)) {
      do.call(rbind,
              lapply(list_draws, function(x) {
                quantile(x$draws_beta_onppseq[, j], probs = quantile_level)
              })
      )
    } else {
      do.call(rbind,
              lapply(list_draws, function(x) {
                quantile(x$draws_sigma_onppseq, probs = quantile_level)
              })
      )
    }
  })
  wis_npp <- unlist(lapply(1:length(q_npp), function(i) {
    mean(wis(rep(true_value[i], N), q_npp[[i]], quantile_level))
  }))
  wis_nppseq <- unlist(lapply(1:length(q_nppseq), function(i) {
    mean(wis(rep(true_value[i], N), q_nppseq[[i]], quantile_level))
  }))
  wis_onpp <- unlist(lapply(1:length(q_onpp), function(i) {
    mean(wis(rep(true_value[i], N), q_onpp[[i]], quantile_level))
  }))
  wis_onppseq <- unlist(lapply(1:length(q_onppseq), function(i) {
    mean(wis(rep(true_value[i], N), q_onppseq[[i]], quantile_level))
  }))
  return(list(hat_wis_npp = wis_npp,
              hat_wis_nppseq = wis_nppseq,
              hat_wis_onpp = wis_onpp,
              hat_wis_onppseq = wis_onppseq
  ))
}

# function to plot BCI
plot_bci <- function(bci, sim, true_val, prob, coverage, scenario, model) {
  bci_npp <- data.frame(sample = 1:length(bci$bci_npp),
                        hattheta = sim$hattheta$hattheta_npp,
                        lower = unlist(lapply(bci$bci_npp, 
                                              function(x) x[[1]])
                        ),
                        upper = unlist(lapply(bci$bci_npp, 
                                              function(x) x[[2]])
                        ),
                        include = unlist(lapply(bci$bci_npp, 
                                                function(x) x[[1]] <= true_val & x[[2]] >= true_val)
                        )
  )
  bci_nppseq <- data.frame(sample = 1:length(bci$bci_nppseq),
                           hattheta = sim$hattheta$hattheta_nppseq,
                           lower = unlist(lapply(bci$bci_nppseq, 
                                                 function(x) x[[1]])
                           ),
                           upper = unlist(lapply(bci$bci_nppseq, 
                                                 function(x) x[[2]])
                           ),
                           include = unlist(lapply(bci$bci_nppseq, 
                                                   function(x) x[[1]] <= true_val & x[[2]] >= true_val)
                           )
  )
  bci_onpp <- data.frame(sample = 1:length(bci$bci_onpp),
                         hattheta = sim$hattheta$hattheta_onpp,
                         lower = unlist(lapply(bci$bci_onpp, 
                                               function(x) x[[1]])
                         ),
                         upper = unlist(lapply(bci$bci_onpp, 
                                               function(x) x[[2]])
                         ),
                         include = unlist(lapply(bci$bci_onpp, 
                                                 function(x) x[[1]] <= true_val & x[[2]] >= true_val)
                         )
  )
  bci_onppseq <- data.frame(sample = 1:length(bci$bci_onppseq),
                            hattheta = sim$hattheta$hattheta_onppseq,
                            lower = unlist(lapply(bci$bci_onppseq, 
                                                  function(x) x[[1]])
                            ),
                            upper = unlist(lapply(bci$bci_onppseq, 
                                                  function(x) x[[2]])
                            ),
                            include = unlist(lapply(bci$bci_onppseq, 
                                                    function(x) x[[1]] <= true_val & x[[2]] >= true_val)
                            )
  )
  
  bci_npp$include <- factor(bci_npp$include, 
                               levels = c(TRUE, FALSE))
  bci_nppseq$include <- factor(bci_nppseq$include, 
                                  levels = c(TRUE, FALSE))
  bci_onpp$include <- factor(bci_onpp$include,
                                levels = c(TRUE, FALSE))
  bci_onppseq$include <- factor(bci_onppseq$include,
                                   levels = c(TRUE, FALSE))
  
  # add dummy row
  bci_npp <- rbind(bci_npp, 
                   data.frame(sample = NA,
                              hattheta = NA,
                              lower = NA,
                              upper = NA,
                              include = factor("FALSE", levels = c("TRUE", "FALSE"))
                   )
  )
  bci_nppseq <- rbind(bci_nppseq, 
                    data.frame(sample = NA,
                               hattheta = NA,
                               lower = NA,
                               upper = NA,
                               include = factor("FALSE", levels = c("TRUE", "FALSE"))
                    )
  )
  bci_onpp <- rbind(bci_onpp, 
                    data.frame(sample = NA,
                               hattheta = NA,
                               lower = NA,
                               upper = NA,
                               include = factor("FALSE", levels = c("TRUE", "FALSE"))
                    )
  )
  bci_onppseq <- rbind(bci_onppseq, 
                       data.frame(sample = NA,
                                  hattheta = NA,
                                  lower = NA,
                                  upper = NA,
                                  include = factor("FALSE", levels = c("TRUE", "FALSE"))
                       )
  )
  
  
  plot_npp <- ggplot() + 
    geom_pointrange(data = bci_npp,
                    mapping = aes(x = sample, y = hattheta,
                                  ymin = lower, ymax = upper, colour = include),
                    na.rm = TRUE,
                    fatten = 0.5) +
    geom_hline(yintercept = 0.5, linetype = "longdash") +
    ggtitle(paste(
      "NPP ", "(",
      formatC(coverage$coverage_npp, 
              format = "f", 
              digits = 2),
      ")",
      sep = ""
    )) + 
    theme_gray() +
    scale_color_manual(values = c("TRUE" = "#00BFC4", "FALSE" = "tomato"), drop = FALSE) +
    ylim(0, 1) +
    xlab("") +
    ylab("")
  
  plot_nppseq <- ggplot() +
    geom_pointrange(data = bci_nppseq,
                    mapping = aes(x = sample, y = hattheta,
                                  ymin = lower, ymax = upper, colour = include),
                    na.rm = TRUE,
                    fatten = 0.5) +
    geom_hline(yintercept = 0.5, linetype = "longdash") +
    ggtitle(paste(
      "NPP-SEQ ", "(",
      formatC(coverage$coverage_nppseq, 
              format = "f", 
              digits = 2), 
      ")",
      sep = ""
    )) +
    theme_gray() +
    scale_color_manual(values = c("TRUE" = "#00BFC4", "FALSE" = "tomato"), drop = FALSE) +
    ylim(0, 1) +
    xlab("") +
    ylab("") +
    scale_y_continuous(limits = c(0, 1), breaks = NULL) 
  
  plot_onpp <- ggplot() +
    geom_pointrange(data = bci_onpp,
                    mapping = aes(x = sample, y = hattheta,
                                  ymin = lower, ymax = upper, colour = include),
                    na.rm = TRUE,
                    fatten = 0.5) +
    geom_hline(yintercept = 0.5, linetype = "longdash") +
    ggtitle(paste(
      "ONPP ", "(",
      formatC(coverage$coverage_onpp, 
              format = "f", 
              digits = 2),
      ")",
      sep = ""
    )) +
    theme_gray() +
    scale_color_manual(values = c("TRUE" = "#00BFC4", "FALSE" = "tomato"), drop = FALSE) +
    ylim(0, 1) +
    xlab("") +
    ylab("") +
    scale_y_continuous(limits = c(0, 1), breaks = NULL) 
  
  plot_onppseq <- ggplot() +
    geom_pointrange(data = bci_onppseq,
                    mapping = aes(x = sample, y = hattheta,
                                  ymin = lower, ymax = upper, colour = include),
                    na.rm = TRUE,
                    fatten = 0.5) +
    geom_hline(yintercept = 0.5, linetype = "longdash") +
    ggtitle(paste(
      "ONPP-SEQ ", "(",
      formatC(coverage$coverage_onppseq, 
              format = "f", 
              digits = 2),
      ")",
      sep = ""
    )) +
    theme_gray() +
    scale_color_manual(values = c("TRUE" = "#00BFC4", "FALSE" = "tomato"), drop = FALSE) +
    ylim(0, 1) +
    xlab("") +
    ylab("") +
    scale_y_continuous(limits = c(0, 1), breaks = NULL) 
  
  # combine plots
  title <- paste(prob*100, "%",
                 " BCI for scenario ", 
                 scenario, 
                 " (", 
                 # model, 
                 ")", 
                 sep = "")
  comb_plot <- (plot_npp + plot_nppseq + plot_onpp + plot_onppseq) + 
    plot_layout(guides = "collect", axis_titles = "collect", ncol = 4) &
    theme(legend.position = "none") &
    theme(text = element_text(size = 16),        
          axis.title = element_text(size = 18),  
          axis.text = element_text(size = 16),   
          legend.title = element_text(size = 18),
          legend.text = element_text(size = 16))
  return(comb_plot)
}

plot_sce_lm <- function(j, hat) {
  sim <- hat[[j]]
  true_value <- c(betastar, sgstar)
  sce <- as.roman(ceiling(j/3))
  sce <- ifelse(j %% 3 == 1, paste0(sce,".I"), ifelse(j %% 3 == 2, paste0(sce,".II"), paste0(sce,".III")))
  plots_delta <- lapply(seq_along(sim$hatdelta), function(i) {
    plot <- plot_boxplot(sim$hatdelta[[i]], bquote(hat(delta)[.(i)])) + coord_flip() +
      ggtitle(paste("Scenario ", sce))
    if (i>1) {
      plot <- plot + scale_x_discrete(labels = NULL) + ggtitle(NULL)
    }
    return(plot)
  })
  plot_theta <- lapply(seq_along(sim$hattheta), function(i) {
    if(i <= length(betastar)) {
      name <- bquote(hat(beta)[.(i)])
      plot <- plot_boxplot(sim$hattheta[[i]], name) + 
        geom_hline(yintercept = true_value[i], linetype = "dotted", color = "red", size = 1) +
        coord_flip()
    } else {
      name <- bquote(hat(sigma^2))
      plot <- plot_boxplot(sim$hattheta[[i]], name) + 
        geom_hline(yintercept = true_value[i]^2, linetype = "dotted", color = "red", size = 1) +
        coord_flip()
    }
    
    if (i>1) {
      plot <- plot + scale_x_discrete(labels = NULL)
    }
    return(plot)
  })
  plot_delta_ <- Reduce(`+`, plots_delta) + plot_layout(ncol = length(plots_delta)) + 
    scale_x_discrete(labels = NULL)
  plot_theta_ <- Reduce(`+`, plot_theta) + plot_layout(ncol = length(plot_theta)) + 
    scale_x_discrete(labels = NULL)
  return(list(plot_delta_,
              plot_theta_))
}


plot_sim_delta <- function(sce, sim, sim_sces){
  name_sce <- as.roman(ceiling(sce/3))
  name_sce <- ifelse(sce %% 3 == 1, paste0(name_sce,".I"), 
                     ifelse(sce %% 3 == 2, paste0(name_sce,".II"), paste0(name_sce,".III")))
  # plot_delta1 <- ggplot() +
  #   geom_density(aes(x = sim_sces[[sce]]$delta[[sim]]$delta_npp[,1], color = "NPP"), linewidth = 0.8) +
  #   geom_density(aes(x = sim_sces[[sce]]$delta[[sim]]$delta_nppseq[,1], color = "NPP_SEQ"), linewidth = 0.8) +
  #   geom_density(aes(x = sim_sces[[sce]]$delta[[sim]]$delta_onpp[,1], color = "ONPP"), linewidth = 0.8) +
  #   geom_density(aes(x = sim_sces[[sce]]$delta[[sim]]$delta_onppseq[,1], color = "ONPP_SEQ"), linewidth = 0.8) +
  #   scale_color_manual(name = NULL, 
  #                      values = c("NPP" = RColorBrewer::brewer.pal(4, "Set1")[1], 
  #                                 "NPP_SEQ" = RColorBrewer::brewer.pal(4, "Set1")[2],
  #                                 "ONPP" = RColorBrewer::brewer.pal(4, "Set1")[3],
  #                                 "ONPP_SEQ" = RColorBrewer::brewer.pal(4, "Set1")[4])) +
  #   labs(x = expression(delta[1]),
  #        y = NULL) +
  #   ggtitle(paste0("Scenario ", name_sce))
  # plot_delta2 <- ggplot() +
  #   geom_density(aes(x = sim_sces[[sce]]$delta[[sim]]$delta_npp[,2], color = "NPP"), linewidth = 0.8) +
  #   geom_density(aes(x = sim_sces[[sce]]$delta[[sim]]$delta_nppseq[,2], color = "NPP_SEQ"), linewidth = 0.8) +
  #   geom_density(aes(x = sim_sces[[sce]]$delta[[sim]]$delta_onpp[,2], color = "ONPP"), linewidth = 0.8) +
  #   geom_density(aes(x = sim_sces[[sce]]$delta[[sim]]$delta_onppseq[,2], color = "ONPP_SEQ"), linewidth = 0.8) +
  #   scale_color_manual(name = NULL, 
  #                      values = c("NPP" = RColorBrewer::brewer.pal(4, "Set1")[1], 
  #                                 "NPP_SEQ" = RColorBrewer::brewer.pal(4, "Set1")[2],
  #                                 "ONPP" = RColorBrewer::brewer.pal(4, "Set1")[3],
  #                                 "ONPP_SEQ" = RColorBrewer::brewer.pal(4, "Set1")[4])) +
  #   labs(x = expression(delta[2]),
  #        y = NULL)
  # plot_delta3 <- ggplot() +
  #   geom_density(aes(x = sim_sces[[sce]]$delta[[sim]]$delta_npp[,3], color = "NPP"), linewidth = 0.8) +
  #   geom_density(aes(x = sim_sces[[sce]]$delta[[sim]]$delta_nppseq[,3], color = "NPP_SEQ"), linewidth = 0.8) +
  #   geom_density(aes(x = sim_sces[[sce]]$delta[[sim]]$delta_onpp[,3], color = "ONPP"), linewidth = 0.8) +
  #   geom_density(aes(x = sim_sces[[sce]]$delta[[sim]]$delta_onppseq[,3], color = "ONPP_SEQ"), linewidth = 0.8) +
  #   scale_color_manual(name = NULL, 
  #                      values = c("NPP" = RColorBrewer::brewer.pal(4, "Set1")[1], 
  #                                 "NPP_SEQ" = RColorBrewer::brewer.pal(4, "Set1")[2],
  #                                 "ONPP" = RColorBrewer::brewer.pal(4, "Set1")[3],
  #                                 "ONPP_SEQ" = RColorBrewer::brewer.pal(4, "Set1")[4])) +
  #   labs(x = expression(delta[3]),
  #        y = NULL)
  # 
  # combined_plots <- (plot_delta1 + plot_delta2 + plot_delta3 + plot_theta) + 
  #   plot_layout(ncol = 4, guides = "collect")
  
  
  plot_npp <- ggplot() +
    geom_density(aes(x = sim_sces[[sce]]$delta[[sim]]$delta_npp[,1], color = "delta1"), linewidth = 0.8) +
    geom_density(aes(x = sim_sces[[sce]]$delta[[sim]]$delta_npp[,2], color = "delta2"), linewidth = 0.8) +
    geom_density(aes(x = sim_sces[[sce]]$delta[[sim]]$delta_npp[,3], color = "delta3"), linewidth = 0.8) +
    scale_color_manual(name = NULL, 
                       values = c("delta1" = RColorBrewer::brewer.pal(3, "Set1")[1], 
                                  "delta2" = RColorBrewer::brewer.pal(3, "Set1")[2],
                                  "delta3" = RColorBrewer::brewer.pal(3, "Set1")[3]),
                       labels = c(expression(delta[1]), 
                                  expression(delta[2]),
                                  expression(delta[3]))) +
    labs(x = expression(delta),
         y = NULL) +
    ggtitle(label = paste0("Scenario ", name_sce),
            subtitle = "NPP")
  
  plot_nppseq <- ggplot() +
    geom_density(aes(x = sim_sces[[sce]]$delta[[sim]]$delta_nppseq[,1], color = "delta1"), linewidth = 0.8) +
    geom_density(aes(x = sim_sces[[sce]]$delta[[sim]]$delta_nppseq[,2], color = "delta2"), linewidth = 0.8) +
    geom_density(aes(x = sim_sces[[sce]]$delta[[sim]]$delta_nppseq[,3], color = "delta3"), linewidth = 0.8) +
    scale_color_manual(name = NULL, 
                       values = c("delta1" = RColorBrewer::brewer.pal(3, "Set1")[1], 
                                  "delta2" = RColorBrewer::brewer.pal(3, "Set1")[2],
                                  "delta3" = RColorBrewer::brewer.pal(3, "Set1")[3]),
                       labels = c(expression(delta[1]), 
                                  expression(delta[2]),
                                  expression(delta[3]))) +
    labs(x = expression(delta),
         y = NULL) +
    ggtitle(label = NULL,
            subtitle =  "NPP-SEQ")
  
  plot_onpp <- ggplot() +
    geom_density(aes(x = sim_sces[[sce]]$delta[[sim]]$delta_onpp[,1], color = "delta1"), linewidth = 0.8) +
    geom_density(aes(x = sim_sces[[sce]]$delta[[sim]]$delta_onpp[,2], color = "delta2"), linewidth = 0.8) +
    geom_density(aes(x = sim_sces[[sce]]$delta[[sim]]$delta_onpp[,3], color = "delta3"), linewidth = 0.8) +
    scale_color_manual(name = NULL, 
                       values = c("delta1" = RColorBrewer::brewer.pal(3, "Set1")[1], 
                                  "delta2" = RColorBrewer::brewer.pal(3, "Set1")[2],
                                  "delta3" = RColorBrewer::brewer.pal(3, "Set1")[3]),
                       labels = c(expression(delta[1]), 
                                  expression(delta[2]),
                                  expression(delta[3]))) +
    labs(x = expression(delta),
         y = NULL) +
    ggtitle(label = NULL,
            subtitle = "ONPP")
  
  plot_onppseq <- ggplot() +
    geom_density(aes(x = sim_sces[[sce]]$delta[[sim]]$delta_onppseq[,1], color = "delta1"), linewidth = 0.8) +
    geom_density(aes(x = sim_sces[[sce]]$delta[[sim]]$delta_onppseq[,2], color = "delta2"), linewidth = 0.8) +
    geom_density(aes(x = sim_sces[[sce]]$delta[[sim]]$delta_onppseq[,3], color = "delta3"), linewidth = 0.8) +
    scale_color_manual(name = NULL, 
                       values = c("delta1" = RColorBrewer::brewer.pal(3, "Set1")[1], 
                                  "delta2" = RColorBrewer::brewer.pal(3, "Set1")[2],
                                  "delta3" = RColorBrewer::brewer.pal(3, "Set1")[3]),
                       labels = c(expression(delta[1]), 
                                  expression(delta[2]),
                                  expression(delta[3]))) +
    labs(x = expression(delta),
         y = NULL) +
    ggtitle(label = NULL,
            subtitle = "ONPP-SEQ")
  combined_plots <- (plot_npp + plot_nppseq + plot_onpp + plot_onppseq) + 
    plot_layout(ncol = 4, guides = "collect")
  return(combined_plots)
}

plot_sim_theta <- function(sce, sim, sim_sces){
  name_sce <- as.roman(ceiling(sce/3))
  name_sce <- ifelse(sce %% 3 == 1, paste0(name_sce,".I"), 
                     ifelse(sce %% 3 == 2, paste0(name_sce,".II"), paste0(name_sce,".III")))
  plot_theta <- ggplot() +
    geom_density(aes(x = sim_sces[[sce]]$theta[[sim]]$theta_npp, color = "NPP"), linewidth = 0.8) +
    geom_density(aes(x = sim_sces[[sce]]$theta[[sim]]$theta_nppseq, color = "NPP_SEQ"), linewidth = 0.8) +
    geom_density(aes(x = sim_sces[[sce]]$theta[[sim]]$theta_onpp, color = "ONPP"), linewidth = 0.8) +
    geom_density(aes(x = sim_sces[[sce]]$theta[[sim]]$theta_onppseq, color = "ONPP_SEQ"), linewidth = 0.8) +
    geom_vline(xintercept = true_value, linetype = "dotted", color = "black", size = 1) +
    scale_color_manual(name = NULL,
                       values = c("NPP" = RColorBrewer::brewer.pal(4, "Set1")[1],
                                  "NPP_SEQ" = RColorBrewer::brewer.pal(4, "Set1")[2],
                                  "ONPP" = RColorBrewer::brewer.pal(4, "Set1")[3],
                                  "ONPP_SEQ" = RColorBrewer::brewer.pal(4, "Set1")[4])) +
    labs(x = expression(theta),
         y = NULL) +
    xlim(0,1) +
    ggtitle(label = paste0("Scenario ", name_sce))
  return(plot_theta)
}

plot_sce3_theta <- function(sub_sce, sim, sim_sces){
  name_sce <- as.roman(ceiling(3))
  name_sce <- ifelse(sub_sce %% 3 == 1, paste0(name_sce,".I"), 
                     ifelse(sub_sce %% 3 == 2, paste0(name_sce,".II"), paste0(name_sce,".III")))
  sim_sub_sce <- sim_sces[[sub_sce]]
  plot_theta <- lapply(1:length(sim_sub_sce), function(i) {
    plot <- ggplot() +
    geom_density(aes(x = sim_sub_sce[[i]]$theta[[sim]]$theta_npp, color = "NPP"), linewidth = 0.8) +
    geom_density(aes(x = sim_sub_sce[[i]]$theta[[sim]]$theta_onpp, color = "ONPP"), linewidth = 0.8) +
    geom_vline(xintercept = true_value, linetype = "dotted", color = "black", size = 1) +
    scale_color_manual(name = NULL,
                       values = c("NPP" = RColorBrewer::brewer.pal(4, "Set2")[1],
                                  "ONPP" = RColorBrewer::brewer.pal(4, "Set2")[2])) +
    labs(x = expression(theta),
         y = NULL) +
    xlim(0,1)
    if (i == 1) {
      plot <- plot +
          ggtitle(label = paste0("Scenario ", name_sce),
                  subtitle = paste0("(", LETTERS[i], ")"))
    } else {
      plot <- plot + 
        ggtitle(label = NULL,
                subtitle = paste0("(", LETTERS[i], ")"))
    }
    return(plot)
  })
  combined_plots <- wrap_plots(plot_theta, ncol = length(plot_theta)) + 
    plot_layout(guides = "collect")
  return(combined_plots)
  }

plot_sce3_eta <- function(sub_sce, sim, sim_sces){
  name_sce <- as.roman(ceiling(3))
  name_sce <- ifelse(sub_sce %% 3 == 1, paste0(name_sce,".I"), 
                     ifelse(sub_sce %% 3 == 2, paste0(name_sce,".II"), paste0(name_sce,".III")))
  
  sim_sub_sce <- sim_sces[[sub_sce]]
  plots_npp <- lapply(1:length(sim_sub_sce), function(i) {
    plot <- ggplot() +
      geom_density(aes(x = sim_sub_sce[[i]]$delta[[sim]]$delta_npp[,1], color = "delta1"), linewidth = 0.8) +
      geom_density(aes(x = sim_sub_sce[[i]]$delta[[sim]]$delta_npp[,2], color = "delta2"), linewidth = 0.8) +
      geom_density(aes(x = sim_sub_sce[[i]]$delta[[sim]]$delta_npp[,3], color = "delta3"), linewidth = 0.8) +
      scale_color_manual(name = NULL, 
                         values = c("delta1" = RColorBrewer::brewer.pal(3, "Set1")[1], 
                                    "delta2" = RColorBrewer::brewer.pal(3, "Set1")[2],
                                    "delta3" = RColorBrewer::brewer.pal(3, "Set1")[3]),
                         labels = c(expression(eta[1]), 
                                    expression(eta[2]),
                                    expression(eta[3]))) +
      labs(x = expression(eta),
           y = NULL)
    if (i == 1) {
      plot <- plot +
          ggtitle(label = paste0("Scenario ", name_sce, "\nNPP"),
                  subtitle = paste0("(", LETTERS[i], ")"))
    } else {
      plot <- plot + 
        ggtitle(label = NULL,
                subtitle = paste0("(", LETTERS[i], ")"))
    }
    return(plot)
    })
  
  plots_onpp <- lapply(1:length(sim_sub_sce), function(i) {
    plot <- ggplot() +
      geom_density(aes(x = sim_sub_sce[[i]]$delta[[sim]]$delta_onpp[,1], color = "delta1"), linewidth = 0.8) +
      geom_density(aes(x = sim_sub_sce[[i]]$delta[[sim]]$delta_onpp[,2], color = "delta2"), linewidth = 0.8) +
      geom_density(aes(x = sim_sub_sce[[i]]$delta[[sim]]$delta_onpp[,3], color = "delta3"), linewidth = 0.8) +
      scale_color_manual(name = NULL, 
                         values = c("delta1" = RColorBrewer::brewer.pal(3, "Set1")[1], 
                                    "delta2" = RColorBrewer::brewer.pal(3, "Set1")[2],
                                    "delta3" = RColorBrewer::brewer.pal(3, "Set1")[3]),
                         labels = c(expression(eta[1]), 
                                    expression(eta[2]),
                                    expression(eta[3]))) +
      labs(x = expression(eta),
           y = NULL)
      if (i == 1) {
        plot <- plot +
          ggtitle(label = "ONPP")
      }
    return(plot)
  })
  
  combined_plots <- wrap_plots(plots_npp, ncol = length(plots_npp)) / 
    wrap_plots(plots_onpp, ncol = length(plots_onpp)) +
    plot_layout(guides = "collect")
  return(combined_plots)
}

plot_eta_prior <- function(samples){
  plots_npp <- lapply(1:length(samples), function(i) {
    plot <- ggplot() +
      geom_density(aes(x = samples[[i]]$delta_npp[,1], color = "delta1"), linewidth = 0.8) +
      geom_density(aes(x = samples[[i]]$delta_npp[,2], color = "delta2"), linewidth = 0.8) +
      geom_density(aes(x = samples[[i]]$delta_npp[,3], color = "delta3"), linewidth = 0.8) +
      scale_color_manual(name = NULL, 
                         values = c("delta1" = RColorBrewer::brewer.pal(3, "Set1")[1], 
                                    "delta2" = RColorBrewer::brewer.pal(3, "Set1")[2],
                                    "delta3" = RColorBrewer::brewer.pal(3, "Set1")[3]),
                         labels = c(expression(eta[1]), 
                                    expression(eta[2]),
                                    expression(eta[3]))) +
      labs(x = expression(eta),
           y = NULL)
    if (i == 1) {
      plot <- plot +
          ggtitle(label = "NPP",
                  subtitle = paste0("(", LETTERS[i], ")"))
    } else {
      plot <- plot + 
        ggtitle(label = NULL,
                subtitle = paste0("(", LETTERS[i], ")"))
    }
    return(plot)
    })
  
  plots_onpp <- lapply(1:length(samples), function(i) {
    plot <- ggplot() +
      geom_density(aes(x = samples[[i]]$delta_onpp[,1], color = "delta1"), linewidth = 0.8) +
      geom_density(aes(x = samples[[i]]$delta_onpp[,2], color = "delta2"), linewidth = 0.8) +
      geom_density(aes(x = samples[[i]]$delta_onpp[,3], color = "delta3"), linewidth = 0.8) +
      scale_color_manual(name = NULL, 
                         values = c("delta1" = RColorBrewer::brewer.pal(3, "Set1")[1], 
                                    "delta2" = RColorBrewer::brewer.pal(3, "Set1")[2],
                                    "delta3" = RColorBrewer::brewer.pal(3, "Set1")[3]),
                         labels = c(expression(eta[1]), 
                                    expression(eta[2]),
                                    expression(eta[3]))) +
      labs(x = expression(eta),
           y = NULL)
      if (i == 1) {
        plot <- plot +
          ggtitle(label = "ONPP")
      }
    return(plot)
  })

  plots_onpp_gamma <- lapply(1:length(samples), function(i) {
    plot <- ggplot() +
      geom_density(aes(x = samples[[i]]$delta_onpp[,1], color = "delta1"), linewidth = 0.8) +
      geom_density(aes(x = samples[[i]]$delta_onpp[,2] - samples[[i]]$delta_onpp[,1], 
        color = "delta2"), linewidth = 0.8) +
      geom_density(aes(x = samples[[i]]$delta_onpp[,3] - samples[[i]]$delta_onpp[,2], 
        color = "delta3"), linewidth = 0.8) +
      scale_color_manual(name = NULL, 
                         values = c("delta1" = RColorBrewer::brewer.pal(3, "Set3")[1], 
                                    "delta2" = RColorBrewer::brewer.pal(3, "Set3")[3],
                                    "delta3" = RColorBrewer::brewer.pal(3, "Set3")[5]),
                         labels = c(expression(gamma[1]), 
                                    expression(gamma[2]),
                                    expression(gamma[3]))) +
      labs(x = expression(gamma),
           y = NULL)
    return(plot)
  })
  
  combined_plots <- wrap_plots(plots_npp, ncol = length(plots_npp)) / 
    wrap_plots(plots_onpp, ncol = length(plots_onpp)) /
    wrap_plots(plots_onpp_gamma, ncol = length(plots_onpp_gamma)) +
    plot_layout(guides = "collect")
  return(combined_plots)
}

plot_theta_prior <- function(samples){
  plot_theta <- lapply(1:length(samples), function(i) {
    plot <- 
    ggplot() +
    geom_density(aes(x = samples[[i]], color = "dens"), linewidth = 0.8) +
    scale_color_manual(name = NULL,
                       values = c("dens" = RColorBrewer::brewer.pal(4, "Set2")[3])) +
    labs(x = expression(theta),
         y = NULL) +
    xlim(0,1) + 
      theme(legend.position = "none") +
        ggtitle(label = NULL,
                subtitle = paste0("(", LETTERS[i], ")"))
    return(plot)
  })
  combined_plots <- wrap_plots(plot_theta, ncol = length(plot_theta)) + 
    plot_layout(guides = "collect")
  return(combined_plots)
  }

