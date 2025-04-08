setwd("C:/Users/Ezequiel/OneDrive - Fundacao Getulio Vargas - FGV/MSC_MAp_CD/onpp-analysis")

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

# run auxiliary functions
source("code/aux_fun_sim.R")

# Compile the model
gamma_model_bin <- cmdstan_model("code/gamma_bin.stan")
delta_model_bin <- cmdstan_model("code/delta_bin.stan")

#-------------------------------------------------------------------------------
# Run simulations for the binomial model with different scenarios
#-------------------------------------------------------------------------------

# define scenarios
sce1 <- list(n0 = c(30, 30, 30),
              n = 30,
              a = 1/2,
              b = 1/2,
              al = 1/2,
              bl = 1/2,
              theta0 = c(0.15, 0.3, 0.45),
              theta = 0.5,
              alpha = rep(1/4, 4)
)
sce2 <- list(n0 = c(30, 30, 30),
             n = 30,
             a = 1/2,
             b = 1/2,
             al = 1/2,
             bl = 1/2,
             theta0 = c(0.2, 0.7, 0.4),
             theta = 0.5,
             alpha = rep(1/4, 4)
)
sce3 <- list(n0 = c(30, 30, 30),
             n = 30,
             a = 1/2,
             b = 1/2,
             al = 1/2,
             bl = 1/2,
             theta0 = c(0.15, 0.85, 0.5),
             theta = 0.5,
             alpha = rep(1/4, 4)
)
sce4 <- list(n0 = c(30, 30, 30),
             n = 30,
             a = 1/2,
             b = 1/2,
             al = 1/2,
             bl = 1/2,
             theta0 = c(0.4, 0.4, 0.4),
             theta = 0.5,
             alpha = rep(1/4, 4)
)
sce5 <- list(n0 = c(30, 30, 30),
             n = 30,
             a = 1/2,
             b = 1/2,
             al = 1/2,
             bl = 1/2,
             theta0 = c(0.1, 0.1, 0.1),
             theta = 0.5,
             alpha = rep(1/4, 4)
)
sce6 <- list(n0 = c(30, 30, 30),
             n = 30,
             a = 1/2,
             b = 1/2,
             al = 1/2,
             bl = 1/2,
             theta0 = c(0.35, 0.2, 0.05),
             theta = 0.5,
             alpha = rep(1/4, 4)
)
sce7 <- list(n0 = c(30, 30, 30),
             n = 30,
             a = 1/2,
             b = 1/2,
             al = 1/2,
             bl = 1/2,
             theta0 = c(0.5, 0.5, 0.5),
             theta = 0.5,
             alpha = rep(1/4, 4)
)

# run simulations
sim1 <- sim_sce(model = "bin",14, 200, sce1, gamma_model_bin, delta_model_bin)
sim2 <- sim_sce(model = "bin",14, 200, sce2, gamma_model_bin, delta_model_bin)
sim3 <- sim_sce(model = "bin",14, 200, sce3, gamma_model_bin, delta_model_bin)
sim4 <- sim_sce(model = "bin",14, 200, sce4, gamma_model_bin, delta_model_bin)
sim5 <- sim_sce(model = "bin",14, 200, sce5, gamma_model_bin, delta_model_bin)
sim6 <- sim_sce(model = "bin",14, 200, sce6, gamma_model_bin, delta_model_bin)
sim7 <- sim_sce(model = "bin",14, 200, sce7, gamma_model_bin, delta_model_bin)

# save results
save(sim1, file = "results/sim1_bin.RData")
save(sim2, file = "results/sim2_bin.RData")
save(sim3, file = "results/sim3_bin.RData")
save(sim4, file = "results/sim4_bin.RData")
save(sim5, file = "results/sim5_bin.RData")
save(sim6, file = "results/sim6_bin.RData")
save(sim7, file = "results/sim7_bin.RData")

#-------------------------------------------------------------------------------
# Run simulations for the binomial model with predefined data
#-------------------------------------------------------------------------------

raw_data1 <- list(n0 = c(10, 10),
             n = 10,
             a = 1/2,
             b = 1/2,
             al = 1/2,
             bl = 1/2,
             theta0 = NULL,
             theta = NULL,
             alpha = rep(1/3, 3),
             z0 = c(5, 5),
             z = 5
)

sample_deltaprior_npp <- delta_model_bin$sample(data = list(a = 1/2, 
                                                        b = 1/2, 
                                                        K = 2, 
                                                        n = 10, 
                                                        z = 5, 
                                                        n0 = c(10,10), 
                                                        z0 = c(5,5), 
                                                        al = 1/2,
                                                        bl = 1/2,
                                                        post = 0,
                                                        seq = 0), 
                                            chains = 4, 
                                            parallel_chains = 4, 
                                            iter_warmup = 2000, 
                                            iter_sampling = 2000,
                                            adapt_delta = 0.98,
                                            refresh = 0
)
sample_deltaprior_onpp <- gamma_model_bin$sample(data = list(a = 1/2, 
                                                         b = 1/2, 
                                                         K = 2, 
                                                         n = 10, 
                                                         z = 5, 
                                                         n0 = c(10,10), 
                                                         z0 = c(5,5), 
                                                         alpha = c(1/(3), 1/(3), 1/(3)),
                                                         post = 0,
                                                         seq = 0), 
                                        chains = 4, 
                                        parallel_chains = 4, 
                                        iter_warmup = 2000, 
                                        iter_sampling = 2000,
                                        adapt_delta = 0.98,
                                        refresh = 0
)

raw_sim1 <- sim_sce(model = "bin", 14, 1, raw_data1, gamma_model_bin, delta_model_bin)
delta_raw_sim1_npp <- read_cmdstan_csv(raw_sim1$output_files[[1]]$npp)$post_warmup_draws %>% 
  as_draws_df() %>% 
  select(contains("delta")) %>% 
  select(-contains("deltal")) %>% 
  mutate(dist = "NPP")
delta_raw_sim1_nppseq <- read_cmdstan_csv(raw_sim1$output_files[[1]]$nppseq)$post_warmup_draws %>% 
  as_draws_df() %>% 
  select(contains("delta")) %>% 
  select(-contains("deltal"))%>% 
  mutate(dist = "NPP-SEQ")
delta_raw_sim1_onpp <- read_cmdstan_csv(raw_sim1$output_files[[1]]$onpp)$post_warmup_draws %>%
  as_draws_df() %>% 
  select(contains("delta")) %>% 
  mutate(dist = "ONPP")
delta_raw_sim1_onppseq <- read_cmdstan_csv(raw_sim1$output_files[[1]]$onppseq)$post_warmup_draws %>%
  as_draws_df() %>% 
  select(contains("delta"))%>% 
  mutate(dist = "ONPP-SEQ")
deltaprior_raw_sim1 <- sample_deltaprior_npp %>% 
  as_draws_df() %>% 
  select(contains("delta")) %>% 
  select(-contains("deltal")) %>% 
  mutate(dist = "Prior")
deltaoprior_raw_sim1 <- sample_deltaprior_onpp %>% 
  as_draws_df() %>% 
  select(contains("delta")) %>% 
  select(-contains("deltal")) %>% 
  mutate(dist = "OPrior")


theta_raw_sim1_npp <- data.frame(theta = raw_sim1$theta[[1]]$theta_npp, dist = "NPP")
theta_raw_sim1_nppseq <- data.frame(theta = raw_sim1$theta[[1]]$theta_nppseq, dist = "NPP-SEQ")
theta_raw_sim1_onpp <- data.frame(theta = raw_sim1$theta[[1]]$theta_onpp, dist = "ONPP")
theta_raw_sim1_onppseq <- data.frame(theta = raw_sim1$theta[[1]]$theta_onppseq, dist = "ONPP-SEQ")
thetaprior_raw_sim1 <- data.frame(theta = rbeta(nrow(theta_raw_sim1_npp), 1/2, 1/2), dist = "Prior")

# Rename columns
colnames(delta_raw_sim1_npp) <- c("delta_1", "delta_2", "dist")
colnames(delta_raw_sim1_nppseq) <- c("delta_1", "delta_2", "dist")
colnames(delta_raw_sim1_onpp) <- c("delta_1", "delta_2", "dist")
colnames(delta_raw_sim1_onppseq) <- c("delta_1", "delta_2", "dist")
colnames(deltaprior_raw_sim1) <- c("delta_1", "delta_2", "dist")
colnames(deltaoprior_raw_sim1) <- c("delta_1", "delta_2", "dist")


# Combine both
df_combined_delta_raw_sim1 <- bind_rows(
  delta_raw_sim1_npp %>% pivot_longer(cols = starts_with("delta"), names_to = "delta_dim", values_to = "value"),
  delta_raw_sim1_nppseq %>% pivot_longer(cols = starts_with("delta"), names_to = "delta_dim", values_to = "value"),
  delta_raw_sim1_onpp %>% pivot_longer(cols = starts_with("delta"), names_to = "delta_dim", values_to = "value"),
  delta_raw_sim1_onppseq %>% pivot_longer(cols = starts_with("delta"), names_to = "delta_dim", values_to = "value"),
  deltaprior_raw_sim1 %>% pivot_longer(cols = starts_with("delta"), names_to = "delta_dim", values_to = "value"),
  deltaoprior_raw_sim1 %>% pivot_longer(cols = starts_with("delta"), names_to = "delta_dim", values_to = "value")
)
df_combined_delta_raw_sim1$delta_dim <- factor(df_combined_delta_raw_sim1$delta_dim, 
                                               labels=c('delta_1'=parse(text=TeX('$\\delta_1$')),
                                                        'delta_2'=parse(text=TeX('$\\delta_2$'))))
df_combined_theta_raw_sim1 <- bind_rows(
  theta_raw_sim1_npp,
  theta_raw_sim1_nppseq,
  theta_raw_sim1_onpp,
  theta_raw_sim1_onppseq,
  thetaprior_raw_sim1
)
colors_delta <- c("OPrior" = "#1F78B4",
                  "Prior" = "#33A02C",
                  "ONPP" = "#E31A1C",
                  "ONPP-SEQ" = "#6A3D9A",
                  "NPP" = "#FF7F00",
                  "NPP-SEQ" = "#A6A6A6"
)

post_delta_raw1_bin <- ggplot(df_combined_delta_raw_sim1, aes(x = value, color = dist)) +
  geom_density(size = 0.8) +
  facet_wrap(~delta_dim, scales = "free", 
             labeller = label_parsed
             ) +  # Separate plots for each dimension
  theme_bw() +
  labs(title = "", x = "", y = "") +
  scale_fill_manual(values = colors_delta) +
  scale_color_manual(values = colors_delta) +
  ylim(0, max(c(density(df_combined_delta_raw_sim1$value)$y, 
                density(df_combined_delta_raw_sim1$value)$y)) * 1.8)

post_theta_raw1_bin <- ggplot(df_combined_theta_raw_sim1, aes(x = theta, color = dist)) +
  geom_density(size = 0.8) +
  theme_bw() +
  labs(title = "", x = expression(lambda), y = "") +
  scale_fill_manual(values = colors_delta) +
  scale_color_manual(values = colors_delta) +
  ylim(0, max(c(density(df_combined_theta_raw_sim1$theta)$y, 
                density(df_combined_theta_raw_sim1$theta)$y)) * 1.3) +
  geom_vline(xintercept = 0.5, linetype = "dotted", color = "black", size = 1)

ggsave("results/figures/post_raw1_bin_delta.png", post_delta_raw1_bin, width = 10, height = 5)
ggsave("results/figures/post_raw1_bin_theta.png", post_theta_raw1_bin, width = 5, height = 5)
