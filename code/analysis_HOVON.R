setwd("C:/Users/Ezequiel/OneDrive - Fundacao Getulio Vargas - FGV/MSC_MAp_CD/onpp-analysis")

# Load libraries
library(cmdstanr)
library(bayesplot)
library(MCMCpack)
library(dplyr)
library(tidyverse)
library(posterior)

# Compile the model
gamma_model_bin <- cmdstan_model("code/gamma_bin.stan")
delta_model_bin <- cmdstan_model("code/delta_bin.stan")

n0 <- c(359, 252, 693, 437)
z0 <- c(279, 208, 598, 358)
n <- 259 
z <- 214
a <- 1/2
b <- 1/2
K <- length(n0)
al <- 1/2
bl <- 1/2
alpha <- rep(1/(K+1), K+1)

data_HOVON_npp <- list(n0 = n0,
                       z0 = z0,
                       n = n,
                       z = z,
                       a = a,
                       b = b,
                       K = K,
                       al = al,
                       bl = bl,
                       post = 1,
                       seq = 0
)
data_HOVON_nppseq <- list(n0 = n0,
                          z0 = z0,
                          n = n,
                          z = z,
                          a = a,
                          b = b,
                          K = K,
                          al = al,
                          bl = bl,
                          post = 1,
                          seq = 1
)
data_HOVON_onpp <- list(n0 = n0,
                        z0 = z0,
                        n = n,
                        z = z,
                        a = a,
                        b = b,
                        K = K,
                        alpha = alpha,
                        post = 1,
                        seq = 0
)
data_HOVON_onppseq <- list(n0 = n0,
                           z0 = z0,
                           n = n,
                           z = z,
                           a = a,
                           b = b,
                           K = K,
                           alpha = alpha,
                           post = 1,
                           seq = 1
)

# Fit the model
sample_delta_npp <- delta_model_bin$sample(data = data_HOVON_npp, 
                                           chains = 4, 
                                           parallel_chains = 4, 
                                           iter_warmup = 2000, 
                                           iter_sampling = 2000
)

sample_delta_nppseq <- delta_model_bin$sample(data = data_HOVON_nppseq, 
                                              chains = 4, 
                                              parallel_chains = 4, 
                                              iter_warmup = 2000, 
                                              iter_sampling = 2000
)

sample_delta_onpp <- gamma_model_bin$sample(data = data_HOVON_onpp, 
                                            chains = 4, 
                                            parallel_chains = 4, 
                                            iter_warmup = 2000, 
                                            iter_sampling = 2000,
                                            adapt_delta = 0.99999
)

sample_delta_onppseq <- gamma_model_bin$sample(data = data_HOVON_onppseq, 
                                               chains = 4, 
                                               parallel_chains = 4, 
                                               iter_warmup = 2000, 
                                               iter_sampling = 2000,
                                               adapt_delta = 0.999
)

# Print summary
sample_delta_npp$summary()
sample_delta_nppseq$summary()
sample_delta_onpp$summary()
sample_delta_onppseq$summary()

# Draws for delta
draws_delta_npp <- sample_delta_npp$draws(variables = "delta") %>% as_draws_matrix()
draws_delta_nppseq <- sample_delta_nppseq$draws(variables = "delta") %>% as_draws_matrix()
draws_delta_onpp <- sample_delta_onpp$draws(variables = "delta") %>% as_draws_matrix()
draws_delta_onppseq <- sample_delta_onppseq$draws(variables = "delta") %>% as_draws_matrix()

delta_1 <- data.frame(draws_delta_npp[,1],
                      draws_delta_nppseq[,1],
                      draws_delta_onpp[,1],
                      draws_delta_onppseq[,1]
)
delta_2 <- data.frame(draws_delta_npp[,2],
                      draws_delta_nppseq[,2],
                      draws_delta_onpp[,2],
                      draws_delta_onppseq[,2]
)
delta_3 <- data.frame(draws_delta_npp[,3],
                      draws_delta_nppseq[,3],
                      draws_delta_onpp[,3],
                      draws_delta_onppseq[,3]
)
delta_4 <- data.frame(draws_delta_npp[,4],
                      draws_delta_nppseq[,4],
                      draws_delta_onpp[,4],
                      draws_delta_onppseq[,4]
)

# Draws for theta
draws_theta_npp <- rbeta(nrow(draws_delta_npp), 
                         z + a + draws_delta_npp %*% z0,
                         n - z + b + draws_delta_npp %*% (n0 - z0)
)
draws_theta_nppseq <- rbeta(nrow(draws_delta_nppseq), 
                            z + K*(a-1) + 1 + draws_delta_nppseq %*% z0,
                            n - z + K*(b-1) + 1 + draws_delta_nppseq %*% (n0 - z0)
)
draws_theta_onpp <- rbeta(nrow(draws_delta_onpp), 
                          z + a + draws_delta_onpp %*% z0,
                          n - z + b + draws_delta_onpp %*% (n0 - z0)
)
draws_theta_onppseq <- rbeta(nrow(draws_delta_onppseq), 
                             z + K*(a-1) + 1 + draws_delta_onppseq %*% z0,
                             n - z + K*(b-1) + 1 + draws_delta_onppseq %*% (n0 - z0)
)

draws_theta <- data.frame(draws_theta_npp,
                          draws_theta_nppseq,
                          draws_theta_onpp,
                          draws_theta_onppseq
)

# Plot boxplot for delta
box_delta_1 <- plot_boxplot(delta_1, expression(delta[1]))
box_delta_2 <- plot_boxplot(delta_2, expression(delta[2]))
box_delta_3 <- plot_boxplot(delta_3, expression(delta[3]))
box_delta_4 <- plot_boxplot(delta_4, expression(delta[4]))

# Combine the plots
box_delta <- (box_delta_1 + box_delta_2) / (box_delta_3 + box_delta_4)

# Plot boxplot for theta
box_theta <- plot_boxplot(draws_theta, expression(theta))

# Save the plots
ggsave("results/figures/box_HOVON_delta.png", box_delta, width = 10, height = 5)
ggsave("results/figures/box_HOVON_theta.png", box_theta, width = 5, height = 5)
