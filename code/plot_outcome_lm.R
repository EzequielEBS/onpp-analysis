# Load libraries
library(ggplot2)
library(viridis)
library(grid)
library(patchwork)

# run auxiliary functions
source("code/aux_fun_sim.R")

# Load data
sim1.1 <- readRDS("results/sim_data/lm/sce_I_I.rds")
sim1.2 <- readRDS("results/sim_data/lm/sce_I_II.rds")
sim1.3 <- readRDS("results/sim_data/lm/sce_I_III.rds")
sim2.1 <- readRDS("results/sim_data/lm/sce_II_I.rds")
sim2.2 <- readRDS("results/sim_data/lm/sce_II_II.rds")
sim2.3 <- readRDS("results/sim_data/lm/sce_II_III.rds")
sim3.1 <- readRDS("results/sim_data/lm/sce_III_I.rds")
sim3.2 <- readRDS("results/sim_data/lm/sce_III_II.rds")
sim3.3 <- readRDS("results/sim_data/lm/sce_III_III.rds")
sim4.1 <- readRDS("results/sim_data/lm/sce_IV_I.rds")
sim4.2 <- readRDS("results/sim_data/lm/sce_IV_II.rds")
sim4.3 <- readRDS("results/sim_data/lm/sce_IV_III.rds")

all_sce_data <- list(sim1.1, sim1.2, sim1.3,
                     sim2.1, sim2.2, sim2.3,
                     sim3.1, sim3.2, sim3.3,
                     sim4.1, sim4.2, sim4.3)

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

plot_mean_dens <- function(i) {
  sce_data <- all_sce_data[[i]]
  c <- cs[[i]]
  hat_y <- lapply(sce_data, function(x) mean(x$y))
  hat_y0 <- lapply(sce_data, function(x) lapply(x$y0, mean))
  hat_y <- unlist(hat_y)
  hat_y0_1 <- unlist(lapply(hat_y0, function(x) x[[1]]))
  hat_y0_2 <- unlist(lapply(hat_y0, function(x) x[[2]]))
  hat_y0_3 <- unlist(lapply(hat_y0, function(x) x[[3]]))
  
  
  sce <- as.roman(ceiling(i/3))
  sce <- ifelse(i %% 3 == 1, paste0(sce,".I"), ifelse(i %% 3 == 2, paste0(sce,".II"), paste0(sce,".III")))
  
  plot <- ggplot() +
    geom_density(aes(x = hat_y0_1, fill = "y0_1"), alpha = 0.5) +
    geom_density(aes(x = hat_y0_2, fill = "y0_2"), alpha = 0.5) +
    geom_density(aes(x = hat_y0_3, fill = "y0_3"), alpha = 0.5) +
    geom_density(aes(x = hat_y, fill = "y"), alpha = 0.5) +
    scale_fill_manual(values = viridis(4),
                      labels = c(expression(hat(y)[01]), 
                                 expression(hat(y)[02]), 
                                 expression(hat(y)[03]),
                                 expression(hat(y))),
                      breaks = c("y0_1", "y0_2", "y0_3", "y"),
                      name = expression(D)) +
    labs(x = "Average outcome", y = "") +
    theme_bw() +
    theme(text = element_text(size = 10),        
          axis.title = element_text(size = 10),  
          axis.text = element_text(size = 10),   
          legend.title = element_text(size = 10),
          legend.text = element_text(size = 10),
          strip.text = element_text(size = 10)) +
    theme(legend.position = c(0.88, 0.75),
          legend.background = element_rect(fill = "white", color = "black")) +
    annotation_custom(
      grid::textGrob(
        label = bquote(c == (.(c[1]) ~","~ 
                           .(c[2]) ~","~
                           .(c[3]))),
        x = unit(0.80, "npc"),  
        y = unit(0.95, "npc"),  
        gp = gpar(
          col = "black",        
          fontsize = 10,
          fill = "white"        
        )
      )
    ) +
    ggtitle(paste0("Scenario ", sce))
  return(plot)
}

plots_average_outcome <- lapply(1:length(all_sce_data), plot_mean_dens)
lapply(1:4, function(i) {
  combined_plot <- patchwork::wrap_plots(plots_average_outcome[((i-1)*3 + 1):((i-1)*3 + 3)], ncol = 3)
  sce <- as.roman(i)
  ggsave(filename = paste0("results/figures/lm/dens_avg_out_", sce, "_lm.pdf"),
         plot = combined_plot,
         width = 10.5, height = 5)
})
