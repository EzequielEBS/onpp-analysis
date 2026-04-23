library(cmdstanr)
library(posterior)
library(dplyr)
library(ggplot2)
library(patchwork)

eta_seq_bin <- cmdstan_model("code/models/eta_seq_bin.stan")
eta_seq_normal <- cmdstan_model("code/models/eta_seq_normal.stan")
gamma_seq_bin <- cmdstan_model("code/models/gamma_seq_bin.stan")

generate_data_bin <- function(sce) {
  n0 <- sce$n0
  theta0 <- sce$theta0
  K <- length(n0)
  z0 <- unlist(lapply(1:K, function(i) rbinom(1, n0[i], theta0[i])))
  list(n0 = n0, theta0 = theta0, K = K, z0 = z0)
}

sample_npp_bin <- function(sce, data) {
  n0 <- data$n0
  theta0 <- data$theta0
  K <- data$K
  z0 <- data$z0
  standata_eta <- list(a = sce$a,
                   b = sce$b,
                   K = K,
                   z0 = z0,
                   n0 = n0,
                   al = sce$al,
                   bl = sce$bl
  )
  sample_eta_nppseq <- eta_seq_bin$sample(data = standata_eta, 
                                          chains = 4, 
                                          parallel_chains = 4, 
                                          iter_warmup = 2000, 
                                          iter_sampling = 2000,
                                          adapt_delta = 0.98,
                                          refresh = 0
  )
  draws_eta_nppseq <- sample_eta_nppseq$draws(variables = "delta") %>% 
    as_draws_matrix()
  draws_eta <- do.call(cbind, lapply(1:K, function(i) rbeta(n = nrow(draws_eta_nppseq), 
                                                          shape1 = sce$al,
                                                          shape2 = sce$bl
  )))
  draws_thetaseq <- rbeta(n = nrow(draws_eta_nppseq), 
                          shape1 = K * (sce$a - 1) + 1,
                          shape2 = K * (sce$b - 1) + 1
  )
  draws_theta <- rbeta(n = nrow(draws_eta_nppseq), 
                       shape1 = sce$a,
                       shape2 = sce$b
  )
  list(
    draws_eta_nppseq = draws_eta_nppseq,
    draws_eta = draws_eta,
    draws_thetaseq = draws_thetaseq,
    draws_theta = draws_theta
  )
}

plot_npp_bin <- function(draws, sce) {
  draws_eta_nppseq <- draws$draws_eta_nppseq
  draws_eta <- draws$draws_eta
  draws_thetaseq <- draws$draws_thetaseq
  draws_theta <- draws$draws_theta
  K <- ncol(draws_eta)

  plot_eta1 <- ggplot() +
    geom_density(aes(x = draws_eta_nppseq[, 1], color = "SEQ"), size = 1) +
    geom_density(aes(x = draws_eta[, 1], color = "pi"), size = 1) +
    labs(x = expression(eta[1]), y = "") +
    scale_color_manual(name = NULL, values = c("SEQ" = "blue", "pi" = "red"), 
                        labels = c("SEQ" = "SEQ",
                                    "pi" = "Default")) +
    theme_minimal() +
    ggtitle(
      bquote(pi[A](eta) %~% prod(Beta(.(sce$al), .(sce$bl)), k==1, K) ~ ","~
             pi[0](theta) %~% Beta(.(sce$a), .(sce$b))
             )
    ) +
    guides(color = "none")
  plot_eta2 <- ggplot() +
    geom_density(aes(x = draws_eta_nppseq[, 2], color = "SEQ"), size = 1) +
    geom_density(aes(x = draws_eta[, 2], color = "pi"), size = 1) +
    labs(x = expression(eta[2]), y = "") +
    scale_color_manual(name = NULL, values = c("SEQ" = "blue", "pi" = "red"), 
                        labels = c("SEQ" = "SEQ",
                                    "pi" = "Default")) +
    theme_minimal() +
    guides(color = "none")
  plot_eta3 <- ggplot() +
    geom_density(aes(x = draws_eta_nppseq[, 3], color = "SEQ"), size = 1) +
    geom_density(aes(x = draws_eta[, 3], color = "pi"), size = 1) +
    labs(x = expression(eta[3]), y = "") +
    scale_color_manual(name = NULL, values = c("SEQ" = "blue", "pi" = "red"), 
                        labels = c("SEQ" = "SEQ",
                                    "pi" = "Default")) +
    theme_minimal() +
    guides(color = "none")
  plot_theta <- ggplot() +
    stat_function(fun = dbeta, args = list(shape1 = sce$a, shape2 = sce$b), 
                  aes(color = "pi"), size = 2) +
    stat_function(fun = dbeta, args = list(shape1 = K * (sce$a - 1) + 1, 
                                            shape2 = K * (sce$b - 1) + 1),
                  aes(color = "SEQ"), size = 1) +
    labs(x = expression(theta), y = "") +
    scale_color_manual(name = NULL, values = c("SEQ" = "blue", "pi" = "red"), 
                        labels = c("SEQ" = "SEQ",
                                    "pi" = "Default")) +
    theme_minimal()
  plot_eta1 + plot_eta2 + plot_eta3 + plot_theta + 
  plot_layout(ncol = 4, guides = "collect") & 
  theme(legend.position = "bottom")
}

# min{a,b} > -1/K + 1
sce1_npp_bin <- list(n0 = c(100, 100, 100),
              a = 1,
              b = 1,
              al = 1/2,
              bl = 1/2,
              theta0 = c(0.1, 0.1, 0.1),
              alpha = c(1/4, 1/4, 1/4, 1/4)
)
sce2_npp_bin <- list(n0 = c(100, 100, 100),
              a = 1,
              b = 1,
              al = 1,
              bl = 1,
              theta0 = c(0.1, 0.1, 0.1),
              alpha = c(1, 1, 1, 1)
)

for (i in 1:10) {
  data_sce1_npp_bin <- generate_data_bin(sce1_npp_bin)
  data_sce2_npp_bin <- generate_data_bin(sce2_npp_bin)
  draws_sce1_npp_bin <- sample_npp_bin(sce1_npp_bin, data_sce1_npp_bin)
  draws_sce2_npp_bin <- sample_npp_bin(sce2_npp_bin, data_sce2_npp_bin)
  plot_sce1_npp_bin <- plot_npp_bin(draws_sce1_npp_bin, sce1_npp_bin)
  plot_sce2_npp_bin <- plot_npp_bin(draws_sce2_npp_bin, sce2_npp_bin)
  print(
  plot_sce1_npp_bin/plot_sce2_npp_bin + plot_layout(guides = "collect") & theme(legend.position = "bottom")
  )
}
ggsave("results/figures/bin/seq_npp_bin.png", width = 14, height = 10)

sample_onpp_bin <- function(sce, data) {
  n0 <- data$n0
  theta0 <- data$theta0
  K <- data$K
  z0 <- data$z0
  alpha <- sce$alpha
  standata_eta_seq <- list(
                  a = sce$a,
                  b = sce$b,
                  alpha = sce$alpha,
                  K = K,
                  z0 = z0,
                  n0 = n0,
                  seq = 1
  )
  standata_eta <- list(
                  a = sce$a,
                  b = sce$b,
                  alpha = sce$alpha,
                  K = K,
                  z0 = z0,
                  n0 = n0,
                  seq = 0
  )

  sample_delta_onppseq <- gamma_seq_bin$sample(data = standata_eta_seq,
                                                chains = 4, 
                                                parallel_chains = 4, 
                                                iter_warmup = 2000, 
                                                iter_sampling = 2000,
                                                adapt_delta = 0.98,
                                                refresh = 0
  )
  sample_delta_onpp <- gamma_seq_bin$sample(data = standata_eta, 
                                                chains = 4, 
                                                parallel_chains = 4, 
                                                iter_warmup = 2000, 
                                                iter_sampling = 2000,
                                                adapt_delta = 0.98,
                                                refresh = 0
  )

  draws_eta_onppseq <- sample_delta_onppseq$draws(variables = "delta") %>% as_draws_matrix()
  draws_eta_onpp <- sample_delta_onpp$draws(variables = "gamma") %>% as_draws_matrix()
  draws_thetaseq <- rbeta(n = nrow(draws_eta_onppseq), 
                          shape1 = K * (sce$a - 1) + 1,
                          shape2 = K * (sce$b - 1) + 1
  )
  draws_theta <- rbeta(n = nrow(draws_eta_onppseq), 
                       shape1 = sce$a,
                       shape2 = sce$b
  )
  list(
    draws_eta_onppseq = draws_eta_onppseq,
    draws_eta_onpp = draws_eta_onpp,
    draws_thetaseq = draws_thetaseq,
    draws_theta = draws_theta
  )
}

plot_onpp_bin <- function(draws, sce) {
  draws_eta_onppseq <- draws$draws_eta_onppseq
  draws_eta_onpp <- draws$draws_eta_onpp
  draws_thetaseq <- draws$draws_thetaseq
  draws_theta <- draws$draws_theta
  K <- ncol(draws_eta_onpp)

  plot_eta1 <- ggplot() +
    geom_density(aes(x = draws_eta_onppseq[, 1], color = "SEQ"), size = 1) +
    geom_density(aes(x = draws_eta_onpp[, 1], color = "pi"), size = 1) +
    labs(x = expression(eta[1]), y = "") +
    scale_color_manual(name = NULL, values = c("SEQ" = "blue", "pi" = "red"), 
                        labels = c("SEQ" = "SEQ",
                                    "pi" = "Default")) +
    theme_minimal() +
    ggtitle(
      bquote(pi[A](gamma) %~% Dirichlet(.(paste(sce$alpha, collapse = ","))) ~ ","~
             pi[0](theta) %~% Beta(.(sce$a), .(sce$b))
             )
    ) +
    guides(color = "none")
  plot_eta2 <- ggplot() +
    geom_density(aes(x = draws_eta_onppseq[, 2], color = "SEQ"), size = 1) +
    geom_density(aes(x = draws_eta_onpp[, 2], color = "pi"), size = 1) +
    labs(x = expression(eta[2]), y = "") +
    scale_color_manual(name = NULL, values = c("SEQ" = "blue", "pi" = "red"), 
                        labels = c("SEQ" = "SEQ",
                                    "pi" = "Default")) +
    theme_minimal() +
    guides(color = "none")
  plot_eta3 <- ggplot() +
    geom_density(aes(x = draws_eta_onppseq[, 3], color = "SEQ"), size = 1) +
    geom_density(aes(x = draws_eta_onpp[, 3], color = "pi"), size = 1) +
    labs(x = expression(eta[3]), y = "") +
    scale_color_manual(name = NULL, values = c("SEQ" = "blue", "pi" = "red"), 
                        labels = c("SEQ" = "SEQ",
                                    "pi" = "Default")) +
    theme_minimal() +
    guides(color = "none")
  plot_theta <- ggplot() +
    stat_function(fun = dbeta, args = list(shape1 = sce$a, shape2 = sce$b), 
                  aes(color = "pi"), size = 2) +
    stat_function(fun = dbeta, args = list(shape1 = K * (sce$a - 1) + 1, 
                                            shape2 = K * (sce$b - 1) + 1),
                  aes(color = "SEQ"), size = 1) +
    labs(x = expression(theta), y = "") +
    scale_color_manual(name = NULL, values = c("SEQ" = "blue", "pi" = "red"), 
                        labels = c("SEQ" = "SEQ",
                                    "pi" = "Default")) +
    theme_minimal()
  plot_eta1 + plot_eta2 + plot_eta3 + plot_theta +
  plot_layout(ncol = 4, guides = "collect") &
  theme(legend.position = "bottom")
}

for (i in 1:10) {
  data_sce1_onpp_bin <- generate_data_bin(sce1_npp_bin)
  data_sce2_onpp_bin <- generate_data_bin(sce2_npp_bin)
  draws_sce1_onpp_bin <- sample_onpp_bin(sce1_npp_bin, data_sce1_onpp_bin)
  draws_sce2_onpp_bin <- sample_onpp_bin(sce2_npp_bin, data_sce2_onpp_bin)
  plot_sce1_onpp_bin <- plot_npp_bin(draws_sce1_onpp_bin, sce1_npp_bin)
  plot_sce2_onpp_bin <- plot_npp_bin(draws_sce2_onpp_bin, sce2_npp_bin)
  print(
  plot_sce1_onpp_bin/plot_sce2_onpp_bin + plot_layout(guides = "collect") & 
    theme(legend.position = "bottom")
  )
}


# -----------------------------------------------------------------------------
# Normal data
# -----------------------------------------------------------------------------

generate_data_normal <- function(sce) {
  n0 <- sce$n0
  theta0 <- sce$theta0
  sigmah <- sce$sigmah
  K <- length(n0)
  y0 <- unlist(lapply(1:K, function(i) rnorm(n0[i], mean = theta0[i], sd = sigmah[i])))
  start_idx <- c(1, cumsum(n0) + 1)[1:K]
  list(n0 = n0, theta0 = theta0, sigmah = sigmah, K = K, y0 = y0, start_idx = start_idx)
}

sample_npp_normal <- function(sce, data) {
  mu0 <- sce$mu0
  sigma0 <- sce$sigma0
  n0 <- data$n0
  theta0 <- data$theta0
  sigmah <- data$sigmah
  K <- data$K
  y0 <- data$y0
  start_idx <- data$start_idx
  standata_eta <- list(
                   mu0 = mu0,
                   sigma0 = sigma0,
                   K = K,
                   y0 = y0,
                   n0 = n0,
                   sigmah = sigmah,
                   start_idx = start_idx,
                   al = sce$al,
                   bl = sce$bl
  )
  sample_eta_nppseq <- eta_seq_normal$sample(data = standata_eta, 
                                          chains = 4, 
                                          parallel_chains = 4, 
                                          iter_warmup = 2000, 
                                          iter_sampling = 2000,
                                          adapt_delta = 0.98,
                                          refresh = 0
  )
  draws_eta_nppseq <- sample_eta_nppseq$draws(variables = "delta") %>% 
    as_draws_matrix()
  draws_eta <- do.call(cbind, lapply(1:K, function(i) rbeta(n = nrow(draws_eta_nppseq), 
                                                          shape1 = sce$al,
                                                          shape2 = sce$bl
  )))
  draws_thetaseq <- rnorm(n = nrow(draws_eta_nppseq), 
                          mean = mu0, 
                          sd = sigma0/sqrt(K)
  )
  draws_theta <- rnorm(n = nrow(draws_eta_nppseq), 
                       mean = mu0, 
                       sd = sigma0
  )
  list(
    draws_eta_nppseq = draws_eta_nppseq,
    draws_eta = draws_eta,
    draws_thetaseq = draws_thetaseq,
    draws_theta = draws_theta
  )
}

plot_npp_normal <- function(draws, sce) {
  draws_eta_nppseq <- draws$draws_eta_nppseq
  draws_eta <- draws$draws_eta
  draws_thetaseq <- draws$draws_thetaseq
  draws_theta <- draws$draws_theta
  K <- ncol(draws_eta)

  plot_eta1 <- ggplot() +
    geom_density(aes(x = draws_eta_nppseq[, 1], color = "SEQ"), size = 1) +
    geom_density(aes(x = draws_eta[, 1], color = "pi"), size = 1) +
    labs(x = expression(eta[1]), y = "") +
    scale_color_manual(name = NULL, values = c("SEQ" = "blue", "pi" = "red"), 
                        labels = c("SEQ" = "SEQ",
                                    "pi" = "Default")) +
    theme_minimal() +
    ggtitle(
      bquote(pi[A](eta) %~% prod(Beta(.(sce$al), .(sce$bl)), k==1, K) ~ ","~
             pi[0](theta) %~% Normal(.(sce$mu0), .(sce$sigma0))
             )
    ) +
    guides(color = "none")
  plot_eta2 <- ggplot() +
    geom_density(aes(x = draws_eta_nppseq[, 2], color = "SEQ"), size = 1) +
    geom_density(aes(x = draws_eta[, 2], color = "pi"), size = 1) +
    labs(x = expression(eta[2]), y = "") +
    scale_color_manual(name = NULL, values = c("SEQ" = "blue", "pi" = "red"), 
                        labels = c("SEQ" = "SEQ",
                                    "pi" = "Default")) +
    theme_minimal() +
    guides(color = "none")
  plot_eta3 <- ggplot() +
    geom_density(aes(x = draws_eta_nppseq[, 3], color = "SEQ"), size = 1) +
    geom_density(aes(x = draws_eta[, 3], color = "pi"), size = 1) +
    labs(x = expression(eta[3]), y = "") +
    scale_color_manual(name = NULL, values = c("SEQ" = "blue", "pi" = "red"), 
                        labels = c("SEQ" = "SEQ",
                                    "pi" = "Default")) +
    theme_minimal() +
    guides(color = "none")
  plot_theta <- ggplot() +
    stat_function(fun = dnorm, args = list(mean = sce$mu0, sd = sce$sigma0), 
                  aes(color = "pi"), size = 2) +
    stat_function(fun = dnorm, args = list(mean = sce$mu0, sd = sce$sigma0/sqrt(K)),
                  aes(color = "SEQ"), size = 1) +
    labs(x = expression(theta), y = "") +
    scale_color_manual(name = NULL, values = c("SEQ" = "blue", "pi" = "red"), 
                        labels = c("SEQ" = "SEQ",
                                    "pi" = "Default")) +
    theme_minimal() +
    xlim(sce$mu0 - 4*sce$sigma0, sce$mu0 + 4*sce$sigma0)
  plot_eta1 + plot_eta2 + plot_eta3 + plot_theta +
  plot_layout(ncol = 4, guides = "collect") & 
  theme(legend.position = "bottom")
}

sce1_npp_normal <- list(n0 = c(100, 100, 100),
              mu0 = 0,
              sigma0 = 1,
              al = 1/2,
              bl = 1/2,
              theta0 = c(-3, -3, 2),
              sigmah = c(1, 1, 1)
)
sce2_npp_normal <- list(n0 = c(100, 100, 100),
              mu0 = 0,
              sigma0 = 1,
              al = 1,
              bl = 1,
              theta0 = c(-3, -3, 2),
              sigmah = c(1, 1, 1)
)

for (i in 1:10) {
  data_sce1_npp_normal <- generate_data_normal(sce1_npp_normal)
  data_sce2_npp_normal <- generate_data_normal(sce2_npp_normal)
  draws_sce1_npp_normal <- sample_npp_normal(sce1_npp_normal, data_sce1_npp_normal)
  draws_sce2_npp_normal <- sample_npp_normal(sce2_npp_normal, data_sce2_npp_normal)
  plot_sce1_npp_normal <- plot_npp_normal(draws_sce1_npp_normal, sce1_npp_normal)
  plot_sce2_npp_normal <- plot_npp_normal(draws_sce2_npp_normal, sce2_npp_normal)
  print(
    plot_sce1_npp_normal/plot_sce2_npp_normal + plot_layout(guides = "collect") & theme(legend.position = "bottom")
  )
}
ggsave("results/figures/normal/seq_npp_normal.png", width = 14, height = 10)
  