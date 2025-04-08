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

# load sample
load("results/sim1_poi.RData")
load("results/sim2_poi.RData")
load("results/sim3_poi.RData")
load("results/sim4_poi.RData")
load("results/sim5_poi.RData")
load("results/sim6_poi.RData")
load("results/sim7_poi.RData")

# Compute MSE
mse1 <- colMeans((sim1$hattheta - 6)^2)
mse2 <- colMeans((sim2$hattheta - 6)^2)
mse3 <- colMeans((sim3$hattheta - 6)^2)
mse4 <- colMeans((sim4$hattheta - 6)^2)
mse5 <- colMeans((sim5$hattheta - 6)^2)
mse6 <- colMeans((sim6$hattheta - 6)^2)
mse7 <- colMeans((sim7$hattheta - 6)^2)

bias1 <- colMeans(sim1$hattheta) - 6
bias2 <- colMeans(sim2$hattheta) - 6
bias3 <- colMeans(sim3$hattheta) - 6
bias4 <- colMeans(sim4$hattheta) - 6
bias5 <- colMeans(sim5$hattheta) - 6
bias6 <- colMeans(sim6$hattheta) - 6
bias7 <- colMeans(sim7$hattheta) - 6

# Compute BCI
bci1_95 <- compute_bci(sim1$theta, 0.95)
bci2_95 <- compute_bci(sim2$theta, 0.95)
bci3_95 <- compute_bci(sim3$theta, 0.95)
bci4_95 <- compute_bci(sim4$theta, 0.95)
bci5_95 <- compute_bci(sim5$theta, 0.95)
bci6_95 <- compute_bci(sim6$theta, 0.95)
bci7_95 <- compute_bci(sim7$theta, 0.95)

bci1_90 <- compute_bci(sim1$theta, 0.90)
bci2_90 <- compute_bci(sim2$theta, 0.90)
bci3_90 <- compute_bci(sim3$theta, 0.90)
bci4_90 <- compute_bci(sim4$theta, 0.90)
bci5_90 <- compute_bci(sim5$theta, 0.90)
bci6_90 <- compute_bci(sim6$theta, 0.90)
bci7_90 <- compute_bci(sim7$theta, 0.90)

bci1_80 <- compute_bci(sim1$theta, 0.80)
bci2_80 <- compute_bci(sim2$theta, 0.80)
bci3_80 <- compute_bci(sim3$theta, 0.80)
bci4_80 <- compute_bci(sim4$theta, 0.80)
bci5_80 <- compute_bci(sim5$theta, 0.80)
bci6_80 <- compute_bci(sim6$theta, 0.80)
bci7_80 <- compute_bci(sim7$theta, 0.80)

bci1_50 <- compute_bci(sim1$theta, 0.50)
bci2_50 <- compute_bci(sim2$theta, 0.50)
bci3_50 <- compute_bci(sim3$theta, 0.50)
bci4_50 <- compute_bci(sim4$theta, 0.50)
bci5_50 <- compute_bci(sim5$theta, 0.50)
bci6_50 <- compute_bci(sim6$theta, 0.50)
bci7_50 <- compute_bci(sim7$theta, 0.50)

# compute coverage
coverage1_95 <- compute_coverage(bci1_95, 6)
coverage2_95 <- compute_coverage(bci2_95, 6)
coverage3_95 <- compute_coverage(bci3_95, 6)
coverage4_95 <- compute_coverage(bci4_95, 6)
coverage5_95 <- compute_coverage(bci5_95, 6)
coverage6_95 <- compute_coverage(bci6_95, 6)
coverage7_95 <- compute_coverage(bci7_95, 6)

coverage1_90 <- compute_coverage(bci1_90, 6)
coverage2_90 <- compute_coverage(bci2_90, 6)
coverage3_90 <- compute_coverage(bci3_90, 6)
coverage4_90 <- compute_coverage(bci4_90, 6)
coverage5_90 <- compute_coverage(bci5_90, 6)
coverage6_90 <- compute_coverage(bci6_90, 6)
coverage7_90 <- compute_coverage(bci7_90, 6)

coverage1_80 <- compute_coverage(bci1_80, 6)
coverage2_80 <- compute_coverage(bci2_80, 6)
coverage3_80 <- compute_coverage(bci3_80, 6)
coverage4_80 <- compute_coverage(bci4_80, 6)
coverage5_80 <- compute_coverage(bci5_80, 6)
coverage6_80 <- compute_coverage(bci6_80, 6)
coverage7_80 <- compute_coverage(bci7_80, 6)

coverage1_50 <- compute_coverage(bci1_50, 6)
coverage2_50 <- compute_coverage(bci2_50, 6)
coverage3_50 <- compute_coverage(bci3_50, 6)
coverage4_50 <- compute_coverage(bci4_50, 6)
coverage5_50 <- compute_coverage(bci5_50, 6)
coverage6_50 <- compute_coverage(bci6_50, 6)
coverage7_50 <- compute_coverage(bci7_50, 6)

# plot BCI
plot_bci1_95 <- plot_bci(bci1_95, sim1, 6, 0.95, coverage1_95, "I", "Poisson")
plot_bci2_95 <- plot_bci(bci2_95, sim2, 6, 0.95, coverage2_95, "II", "Poisson")
plot_bci3_95 <- plot_bci(bci3_95, sim3, 6, 0.95, coverage3_95, "III", "Poisson")
plot_bci4_95 <- plot_bci(bci4_95, sim4, 6, 0.95, coverage4_95, "IV", "Poisson")
plot_bci5_95 <- plot_bci(bci5_95, sim5, 6, 0.95, coverage5_95, "V", "Poisson")
plot_bci6_95 <- plot_bci(bci6_95, sim6, 6, 0.95, coverage6_95, "VI", "Poisson")
plot_bci7_95 <- plot_bci(bci7_95, sim7, 6, 0.95, coverage7_95, "VII", "Poisson")

plot_bci1_90 <- plot_bci(bci1_90, sim1, 6, 0.90, coverage1_90, "I", "Poisson")
plot_bci2_90 <- plot_bci(bci2_90, sim2, 6, 0.90, coverage2_90, "II", "Poisson")
plot_bci3_90 <- plot_bci(bci3_90, sim3, 6, 0.90, coverage3_90, "III", "Poisson")
plot_bci4_90 <- plot_bci(bci4_90, sim4, 6, 0.90, coverage4_90, "IV", "Poisson")
plot_bci5_90 <- plot_bci(bci5_90, sim5, 6, 0.90, coverage5_90, "V", "Poisson")
plot_bci6_90 <- plot_bci(bci6_90, sim6, 6, 0.90, coverage6_90, "VI", "Poisson")
plot_bci7_90 <- plot_bci(bci7_90, sim7, 6, 0.90, coverage7_90, "VII", "Poisson")

plot_bci1_80 <- plot_bci(bci1_80, sim1, 6, 0.80, coverage1_80, "I", "Poisson")
plot_bci2_80 <- plot_bci(bci2_80, sim2, 6, 0.80, coverage2_80, "II", "Poisson")
plot_bci3_80 <- plot_bci(bci3_80, sim3, 6, 0.80, coverage3_80, "III", "Poisson")
plot_bci4_80 <- plot_bci(bci4_80, sim4, 6, 0.80, coverage4_80, "IV", "Poisson")
plot_bci5_80 <- plot_bci(bci5_80, sim5, 6, 0.80, coverage5_80, "V", "Poisson")
plot_bci6_80 <- plot_bci(bci6_80, sim6, 6, 0.80, coverage6_80, "VI", "Poisson")
plot_bci7_80 <- plot_bci(bci7_80, sim7, 6, 0.80, coverage7_80, "VII", "Poisson")

plot_bci1_50 <- plot_bci(bci1_50, sim1, 6, 0.50, coverage1_50, "I", "Poisson")
plot_bci2_50 <- plot_bci(bci2_50, sim2, 6, 0.50, coverage2_50, "II", "Poisson")
plot_bci3_50 <- plot_bci(bci3_50, sim3, 6, 0.50, coverage3_50, "III", "Poisson")
plot_bci4_50 <- plot_bci(bci4_50, sim4, 6, 0.50, coverage4_50, "IV", "Poisson")
plot_bci5_50 <- plot_bci(bci5_50, sim5, 6, 0.50, coverage5_50, "V", "Poisson")
plot_bci6_50 <- plot_bci(bci6_50, sim6, 6, 0.50, coverage6_50, "VI", "Poisson")
plot_bci7_50 <- plot_bci(bci7_50, sim7, 6, 0.50, coverage7_50, "VII", "Poisson")

# Save the plots
ggsave("results/figures/bci1_95_poi.png", plot_bci1_95, width = 20, height = 10, dpi = 300)
ggsave("results/figures/bci2_95_poi.png", plot_bci2_95, width = 20, height = 10, dpi = 300)
ggsave("results/figures/bci3_95_poi.png", plot_bci3_95, width = 20, height = 10, dpi = 300)
ggsave("results/figures/bci4_95_poi.png", plot_bci4_95, width = 20, height = 10, dpi = 300)
ggsave("results/figures/bci5_95_poi.png", plot_bci5_95, width = 20, height = 10, dpi = 300)
ggsave("results/figures/bci6_95_poi.png", plot_bci6_95, width = 20, height = 10, dpi = 300)
ggsave("results/figures/bci7_95_poi.png", plot_bci7_95, width = 20, height = 10, dpi = 300)

ggsave("results/figures/bci1_90_poi.png", plot_bci1_90, width = 20, height = 10, dpi = 300)
ggsave("results/figures/bci2_90_poi.png", plot_bci2_90, width = 20, height = 10, dpi = 300)
ggsave("results/figures/bci3_90_poi.png", plot_bci3_90, width = 20, height = 10, dpi = 300)
ggsave("results/figures/bci4_90_poi.png", plot_bci4_90, width = 20, height = 10, dpi = 300)
ggsave("results/figures/bci5_90_poi.png", plot_bci5_90, width = 20, height = 10, dpi = 300)
ggsave("results/figures/bci6_90_poi.png", plot_bci6_90, width = 20, height = 10, dpi = 300)
ggsave("results/figures/bci7_90_poi.png", plot_bci7_90, width = 20, height = 10, dpi = 300)

ggsave("results/figures/bci1_80_poi.png", plot_bci1_80, width = 20, height = 10, dpi = 300)
ggsave("results/figures/bci2_80_poi.png", plot_bci2_80, width = 20, height = 10, dpi = 300)
ggsave("results/figures/bci3_80_poi.png", plot_bci3_80, width = 20, height = 10, dpi = 300)
ggsave("results/figures/bci4_80_poi.png", plot_bci4_80, width = 20, height = 10, dpi = 300)
ggsave("results/figures/bci5_80_poi.png", plot_bci5_80, width = 20, height = 10, dpi = 300)
ggsave("results/figures/bci6_80_poi.png", plot_bci6_80, width = 20, height = 10, dpi = 300)
ggsave("results/figures/bci7_80_poi.png", plot_bci7_80, width = 20, height = 10, dpi = 300)

ggsave("results/figures/bci1_50_poi.png", plot_bci1_50, width = 20, height = 10, dpi = 300)
ggsave("results/figures/bci2_50_poi.png", plot_bci2_50, width = 20, height = 10, dpi = 300)
ggsave("results/figures/bci3_50_poi.png", plot_bci3_50, width = 20, height = 10, dpi = 300)
ggsave("results/figures/bci4_50_poi.png", plot_bci4_50, width = 20, height = 10, dpi = 300)
ggsave("results/figures/bci5_50_poi.png", plot_bci5_50, width = 20, height = 10, dpi = 300)
ggsave("results/figures/bci6_50_poi.png", plot_bci6_50, width = 20, height = 10, dpi = 300)
ggsave("results/figures/bci7_50_poi.png", plot_bci7_50, width = 20, height = 10, dpi = 300)


# ------------------------------------------------------------------------------
# Summary results
# ------------------------------------------------------------------------------


table1 <- data.frame(scenario = rep("I", 4),
                     model = c("NPP", "NPP-SEQ", "ONPP", "ONPP-SEQ"),
                     bias = round(unlist(bias1), 3),
                     mse = round(unlist(mse1), 3),
                     coverage_95 = unlist(coverage1_95),
                     coverage_90 = unlist(coverage1_90),
                     coverage_80 = unlist(coverage1_80),
                     coverage_50 = unlist(coverage1_50)
)
table2 <- data.frame(scenario = rep("II", 4),
                     model = c("NPP", "NPP-SEQ", "ONPP", "ONPP-SEQ"),
                     bias = round(unlist(bias2), 3),
                     mse = round(unlist(mse2), 3),
                     coverage_95 = unlist(coverage2_95),
                     coverage_90 = unlist(coverage2_90),
                     coverage_80 = unlist(coverage2_80),
                     coverage_50 = unlist(coverage2_50)
)
table3 <- data.frame(scenario = rep("III", 4),
                     model = c("NPP", "NPP-SEQ", "ONPP", "ONPP-SEQ"),
                     bias = round(unlist(bias3), 3),
                     mse = round(unlist(mse3), 3),
                     coverage_95 = unlist(coverage3_95),
                     coverage_90 = unlist(coverage3_90),
                     coverage_80 = unlist(coverage3_80),
                     coverage_50 = unlist(coverage3_50)
)
table4 <- data.frame(scenario = rep("IV", 4),
                     model = c("NPP", "NPP-SEQ", "ONPP", "ONPP-SEQ"),
                     bias = round(unlist(bias4), 3),
                     mse = round(unlist(mse4), 3),
                     coverage_95 = unlist(coverage4_95),
                     coverage_90 = unlist(coverage4_90),
                     coverage_80 = unlist(coverage4_80),
                     coverage_50 = unlist(coverage4_50)
)
table5 <- data.frame(scenario = rep("V", 4),
                     model = c("NPP", "NPP-SEQ", "ONPP", "ONPP-SEQ"),
                     bias = round(unlist(bias5), 3),
                     mse = round(unlist(mse5), 3),
                     coverage_95 = unlist(coverage5_95),
                     coverage_90 = unlist(coverage5_90),
                     coverage_80 = unlist(coverage5_80),
                     coverage_50 = unlist(coverage5_50)
)
table6 <- data.frame(scenario = rep("VI", 4),
                     model = c("NPP", "NPP-SEQ", "ONPP", "ONPP-SEQ"),
                     bias = round(unlist(bias6), 3),
                     mse = round(unlist(mse6), 3),
                     coverage_95 = unlist(coverage6_95),
                     coverage_90 = unlist(coverage6_90),
                     coverage_80 = unlist(coverage6_80),
                     coverage_50 = unlist(coverage6_50)
)
table7 <- data.frame(scenario = rep("VII", 4),
                     model = c("NPP", "NPP-SEQ", "ONPP", "ONPP-SEQ"),
                     bias = round(unlist(bias7), 3),
                     mse = round(unlist(mse7), 3),
                     coverage_95 = unlist(coverage7_95),
                     coverage_90 = unlist(coverage7_90),
                     coverage_80 = unlist(coverage7_80),
                     coverage_50 = unlist(coverage7_50)
)

summary_table <- bind_rows(table1, table2, table3, table4, table5, table6, table7)
colnames(summary_table) <- c("Scenario", "Prior", "Bias", "MSE", 
                             "Cov_BCI95", "Cov_BCI90", "Cov_BCI80", "Cov_BCI50")
row.names(summary_table) <- NULL

print(xtable::xtable(summary_table %>% 
                       mutate(
                         MSE = formatC(MSE, format = "f", digits = 3),
                         Bias = formatC(Bias, format = "f", digits = 3)
                       )
), 
type = "latex", include.rownames = FALSE)

# ------------------------------------------------------------------------------
# Plot boxplots
# ------------------------------------------------------------------------------

# Boxplots for delta
box_hatdelta1_1 <- plot_boxplot(sim1$hatdelta[[1]], expression(hat(delta)[1])) + 
  ylim(0, 1)
box_hatdelta1_2 <- plot_boxplot(sim1$hatdelta[[2]], expression(hat(delta)[2])) + 
  ylim(0, 1)
box_hatdelta1_3 <- plot_boxplot(sim1$hatdelta[[3]], expression(hat(delta)[3])) + 
  ylim(0, 1)

box_hatdelta2_1 <- plot_boxplot(sim2$hatdelta[[1]], expression(hat(delta)[1])) + 
  ylim(0, 1)
box_hatdelta2_2 <- plot_boxplot(sim2$hatdelta[[2]], expression(hat(delta)[2])) + 
  ylim(0, 1)
box_hatdelta2_3 <- plot_boxplot(sim2$hatdelta[[3]], expression(hat(delta)[3])) + 
  ylim(0, 1)

box_hatdelta3_1 <- plot_boxplot(sim3$hatdelta[[1]], expression(hat(delta)[1])) + 
  ylim(0, 1)
box_hatdelta3_2 <- plot_boxplot(sim3$hatdelta[[2]], expression(hat(delta)[2])) + 
  ylim(0, 1)
box_hatdelta3_3 <- plot_boxplot(sim3$hatdelta[[3]], expression(hat(delta)[3])) + 
  ylim(0, 1)

box_hatdelta4_1 <- plot_boxplot(sim4$hatdelta[[1]], expression(hat(delta)[1])) + 
  ylim(0, 1)
box_hatdelta4_2 <- plot_boxplot(sim4$hatdelta[[2]], expression(hat(delta)[2])) + 
  ylim(0, 1)
box_hatdelta4_3 <- plot_boxplot(sim4$hatdelta[[3]], expression(hat(delta)[3])) + 
  ylim(0, 1)

box_hatdelta5_1 <- plot_boxplot(sim5$hatdelta[[1]], expression(hat(delta)[1])) + 
  ylim(0, 1)
box_hatdelta5_2 <- plot_boxplot(sim5$hatdelta[[2]], expression(hat(delta)[2])) + 
  ylim(0, 1)
box_hatdelta5_3 <- plot_boxplot(sim5$hatdelta[[3]], expression(hat(delta)[3])) + 
  ylim(0, 1)

box_hatdelta6_1 <- plot_boxplot(sim6$hatdelta[[1]], expression(hat(delta)[1])) + 
  ylim(0, 1)
box_hatdelta6_2 <- plot_boxplot(sim6$hatdelta[[2]], expression(hat(delta)[2])) + 
  ylim(0, 1)
box_hatdelta6_3 <- plot_boxplot(sim6$hatdelta[[3]], expression(hat(delta)[3])) + 
  ylim(0, 1)

box_hatdelta7_1 <- plot_boxplot(sim7$hatdelta[[1]], expression(hat(delta)[1])) + 
  ylim(0, 1)
box_hatdelta7_2 <- plot_boxplot(sim7$hatdelta[[2]], expression(hat(delta)[2])) + 
  ylim(0, 1)
box_hatdelta7_3 <- plot_boxplot(sim7$hatdelta[[3]], expression(hat(delta)[3])) + 
  ylim(0, 1)

# Make a grid of boxplots
box_hatdelta1 <- box_hatdelta1_1 + box_hatdelta1_2 + box_hatdelta1_3
box_hatdelta2 <- box_hatdelta2_1 + box_hatdelta2_2 + box_hatdelta2_3
box_hatdelta3 <- box_hatdelta3_1 + box_hatdelta3_2 + box_hatdelta3_3
box_hatdelta4 <- box_hatdelta4_1 + box_hatdelta4_2 + box_hatdelta4_3
box_hatdelta5 <- box_hatdelta5_1 + box_hatdelta5_2 + box_hatdelta5_3
box_hatdelta6 <- box_hatdelta6_1 + box_hatdelta6_2 + box_hatdelta6_3
box_hatdelta7 <- box_hatdelta7_1 + box_hatdelta7_2 + box_hatdelta7_3

# Boxplots for theta
box_hattheta1 <- plot_boxplot(sim1$hattheta, expression(hat(lambda))) + 
  geom_hline(yintercept = 6, linetype = "dotted", color = "black", size = 1)
box_hattheta2 <- plot_boxplot(sim2$hattheta, expression(hat(lambda))) + 
  geom_hline(yintercept = 6, linetype = "dotted", color = "black", size = 1)
box_hattheta3 <- plot_boxplot(sim3$hattheta, expression(hat(lambda))) +
  geom_hline(yintercept = 6, linetype = "dotted", color = "black", size = 1)
box_hattheta4 <- plot_boxplot(sim4$hattheta, expression(hat(lambda))) +
  geom_hline(yintercept = 6, linetype = "dotted", color = "black", size = 1)
box_hattheta5 <- plot_boxplot(sim5$hattheta, expression(hat(lambda))) +
  geom_hline(yintercept = 6, linetype = "dotted", color = "black", size = 1)
box_hattheta6 <- plot_boxplot(sim6$hattheta, expression(hat(lambda))) +
  geom_hline(yintercept = 6, linetype = "dotted", color = "black", size = 1)
box_hattheta7 <- plot_boxplot(sim7$hattheta, expression(hat(lambda))) +
  geom_hline(yintercept = 6, linetype = "dotted", color = "black", size = 1)

# save the plots
ggsave("results/figures/box_simpoi_hatdelta1.png", box_hatdelta1, width = 10, height = 5)
ggsave("results/figures/box_simpoi_hatdelta2.png", box_hatdelta2, width = 10, height = 5)
ggsave("results/figures/box_simpoi_hatdelta3.png", box_hatdelta3, width = 10, height = 5)
ggsave("results/figures/box_simpoi_hatdelta4.png", box_hatdelta4, width = 10, height = 5)
ggsave("results/figures/box_simpoi_hatdelta5.png", box_hatdelta5, width = 10, height = 5)
ggsave("results/figures/box_simpoi_hatdelta6.png", box_hatdelta6, width = 10, height = 5)
ggsave("results/figures/box_simpoi_hatdelta7.png", box_hatdelta7, width = 10, height = 5)

ggsave("results/figures/box_simpoi_hattheta1.png", box_hattheta1, width = 5, height = 5)
ggsave("results/figures/box_simpoi_hattheta2.png", box_hattheta2, width = 5, height = 5)
ggsave("results/figures/box_simpoi_hattheta3.png", box_hattheta3, width = 5, height = 5)
ggsave("results/figures/box_simpoi_hattheta4.png", box_hattheta4, width = 5, height = 5)
ggsave("results/figures/box_simpoi_hattheta5.png", box_hattheta5, width = 5, height = 5)
ggsave("results/figures/box_simpoi_hattheta6.png", box_hattheta6, width = 5, height = 5)
ggsave("results/figures/box_simpoi_hattheta7.png", box_hattheta7, width = 5, height = 5)