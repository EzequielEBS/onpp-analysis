

source("code/aux_fun_sim.R")
source("code/sim_exp/bin/scenarios_bin.R")

# load sample
sim1.1 <- qs2::qs_read("results/samples/bin/sim1_1_bin.qs2")
sim1.2 <- qs2::qs_read("results/samples/bin/sim1_2_bin.qs2")
sim1.3 <- qs2::qs_read("results/samples/bin/sim1_3_bin.qs2")
sim2.1 <- qs2::qs_read("results/samples/bin/sim2_1_bin.qs2")
sim2.2 <- qs2::qs_read("results/samples/bin/sim2_2_bin.qs2")
sim2.3 <- qs2::qs_read("results/samples/bin/sim2_3_bin.qs2")
sim3.1 <- qs2::qs_read("results/samples/bin/sim3_1_bin.qs2")
sim3.2 <- qs2::qs_read("results/samples/bin/sim3_2_bin.qs2")
sim3.3 <- qs2::qs_read("results/samples/bin/sim3_3_bin.qs2")
sim4.1 <- qs2::qs_read("results/samples/bin/sim4_1_bin.qs2")
sim4.2 <- qs2::qs_read("results/samples/bin/sim4_2_bin.qs2")
sim4.3 <- qs2::qs_read("results/samples/bin/sim4_3_bin.qs2")
sim5.1 <- qs2::qs_read("results/samples/bin/sim5_1_bin.qs2")
sim5.2 <- qs2::qs_read("results/samples/bin/sim5_2_bin.qs2")

sim3.1_eta_flat <- qs2::qs_read("results/samples/bin/sim3_1_bin_eta_flat.qs2")
sim3.2_eta_flat <- qs2::qs_read("results/samples/bin/sim3_2_bin_eta_flat.qs2")
sim3.3_eta_flat <- qs2::qs_read("results/samples/bin/sim3_3_bin_eta_flat.qs2")

sim3.1_eta_gamma_flat <- qs2::qs_read("results/samples/bin/sim3_1_bin_eta_gamma_flat.qs2")
sim3.2_eta_gamma_flat <- qs2::qs_read("results/samples/bin/sim3_2_bin_eta_gamma_flat.qs2")
sim3.3_eta_gamma_flat <- qs2::qs_read("results/samples/bin/sim3_3_bin_eta_gamma_flat.qs2")

sim3.1_all_flat <- qs2::qs_read("results/samples/bin/sim3_1_bin_all_flat.qs2")
sim3.2_all_flat <- qs2::qs_read("results/samples/bin/sim3_2_bin_all_flat.qs2")
sim3.3_all_flat <- qs2::qs_read("results/samples/bin/sim3_3_bin_all_flat.qs2")

samples_prior_sce3.1_all_flat <- qs2::qs_read("results/samples/bin/samples_prior_sce3_1_all_flat.qs2")
samples_prior_sce3.2_all_flat <- qs2::qs_read("results/samples/bin/samples_prior_sce3_2_all_flat.qs2")
samples_prior_sce3.3_all_flat <- qs2::qs_read("results/samples/bin/samples_prior_sce3_3_all_flat.qs2")

# true value
true_value <- 0.5

# combine simulations
sim_sces <- list(
  sim1.1, sim1.2, sim1.3,
  sim2.1, sim2.2, sim2.3,
  sim3.1, sim3.2, sim3.3,
  sim4.1, sim4.2, sim4.3,
  sim5.1, sim5.2
)

plot_eta_theta_density(sim_sces[[14]], replicate = 1, family = "bin",
                        true_theta = 0.5,
                        scenario_label = "Scenario V.I")

plot_eta_theta_boxplot(sim_sces[[13]], family = "bin", true_theta = 0.5,
                        scenario_label = "Scenario V.I")

plot_eta_theta_boxplot(sim_sces_beta_n0_ge_n[[2]], family = "lm",
                        true_theta = c(betastar, sgstar))