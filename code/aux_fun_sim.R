# function to simulate data and run the models
sample_sce_bin <- function(par_list, gamma_model, delta_model){
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
                    post = 1,
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
                       post = 1,
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
                   post = 1,
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
                      post = 1,
                      seq = 1
  )
  # sample from the models
  sample_delta_onpp <- gamma_model$sample(data = data_onpp, 
                                              chains = 4, 
                                              parallel_chains = 4, 
                                              iter_warmup = 2000, 
                                              iter_sampling = 2000,
                                              adapt_delta = 0.98,
                                              refresh = 0
  )
  sample_delta_onppseq <- gamma_model$sample(data = data_onppseq, 
                                                 chains = 4, 
                                                 parallel_chains = 4, 
                                                 iter_warmup = 2000, 
                                                 iter_sampling = 2000,
                                                 adapt_delta = 0.98,
                                                 refresh = 0
  )
  sample_delta_npp <- delta_model$sample(data = data_npp, 
                                             chains = 4, 
                                             parallel_chains = 4, 
                                             iter_warmup = 2000, 
                                             iter_sampling = 2000,
                                             adapt_delta = 0.98,
                                             refresh = 0
  )
  sample_delta_nppseq <- delta_model$sample(data = data_nppseq, 
                                                chains = 4, 
                                                parallel_chains = 4, 
                                                iter_warmup = 2000, 
                                                iter_sampling = 2000,
                                                adapt_delta = 0.98,
                                                refresh = 0
  )
  # save output files
  save_dir <- "results/cmdstan_outputs"
  if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
  sample_delta_npp$save_output_files(dir = save_dir)
  sample_delta_nppseq$save_output_files(dir = save_dir)
  sample_delta_onpp$save_output_files(dir = save_dir)
  sample_delta_onppseq$save_output_files(dir = save_dir)
  
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
              output_files = list(
                npp = sample_delta_npp$output_files(),
                nppseq = sample_delta_nppseq$output_files(),
                onpp = sample_delta_onpp$output_files(),
                onppseq = sample_delta_onppseq$output_files()
              ),
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
                                          parallel_chains = 4, 
                                          iter_warmup = 2000, 
                                          iter_sampling = 2000,
                                          adapt_delta = 0.98,
                                          refresh = 0
  )
  sample_delta_onppseq <- gamma_model$sample(data = data_onppseq, 
                                             chains = 4, 
                                             parallel_chains = 4, 
                                             iter_warmup = 2000, 
                                             iter_sampling = 2000,
                                             adapt_delta = 0.98,
                                             refresh = 0
  )
  sample_delta_npp <- delta_model$sample(data = data_npp, 
                                         chains = 4, 
                                         parallel_chains = 4, 
                                         iter_warmup = 2000, 
                                         iter_sampling = 2000,
                                         adapt_delta = 0.98,
                                         refresh = 0
  )
  sample_delta_nppseq <- delta_model$sample(data = data_nppseq, 
                                            chains = 4, 
                                            parallel_chains = 4, 
                                            iter_warmup = 2000, 
                                            iter_sampling = 2000,
                                            adapt_delta = 0.98,
                                            refresh = 0
  )
  # save output files
  save_dir <- "results/cmdstan_outputs"
  if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
  sample_delta_npp$save_output_files(dir = save_dir)
  sample_delta_nppseq$save_output_files(dir = save_dir)
  sample_delta_onpp$save_output_files(dir = save_dir)
  sample_delta_onppseq$save_output_files(dir = save_dir)
  
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
              output_files = list(
                npp = sample_delta_npp$output_files(),
                nppseq = sample_delta_nppseq$output_files(),
                onpp = sample_delta_onpp$output_files(),
                onppseq = sample_delta_onppseq$output_files()
              ),
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
    library(bayesplot)
    library(MCMCpack)
    library(dplyr)
    library(tidyverse)
    library(posterior)
    library(patchwork)
  })
  # Export the models to the cluster
  clusterExport(cl, 
                varlist = c("gamma_model", 
                            "delta_model"),
                envir = environment()
  )
  
  # Generate a list of parameters for each simulation
  rep_sce <- replicate(num_sim, sce, simplify = FALSE)
  
  # Run simulations in parallel
  if (model == "bin"){
    sample_sce <- sample_sce_bin
  } else if (model == "poi"){
    sample_sce <- sample_sce_poi
  }
  results <- parLapply(cl, 
                       rep_sce, 
                       sample_sce,
                       gamma_model = gamma_model,
                       delta_model = delta_model)
  
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
              output_files = lapply(results, function(x) x$output_files),
              theta = lapply(results, function(x) list(theta_npp = x$theta_npp,
                                                        theta_nppseq = x$theta_nppseq,
                                                        theta_onpp = x$theta_onpp,
                                                        theta_onppseq = x$theta_onppseq)
              )
  ))
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
    geom_boxplot(outlier.colour="black", outlier.shape= NA, 
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
                    na.rm = TRUE) +
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
                    na.rm = TRUE) +
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
    ylab("")
  
  plot_onpp <- ggplot() +
    geom_pointrange(data = bci_onpp,
                    mapping = aes(x = sample, y = hattheta,
                                  ymin = lower, ymax = upper, colour = include),
                    na.rm = TRUE) +
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
    ylab("")
  
  plot_onppseq <- ggplot() +
    geom_pointrange(data = bci_onppseq,
                    mapping = aes(x = sample, y = hattheta,
                                  ymin = lower, ymax = upper, colour = include),
                    na.rm = TRUE) +
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
    ylab("")
  
  # combine plots
  title <- paste(prob*100, "%",
                 " BCI for scenario ", 
                 scenario, 
                 " (", 
                 model, 
                 ")", 
                 sep = "")
  comb_plot <- (plot_npp + plot_nppseq) / (plot_onpp + plot_onppseq) + 
    plot_layout(guides = "collect", axis_titles = "collect") +
    plot_annotation(title = title)
  return(comb_plot)
}
