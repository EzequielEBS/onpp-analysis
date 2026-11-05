# =============================================================================
# Schematic figure: how each Order category (Compatible / Semi-compat. /
# Incompatible / Neutral / Cancellation -- the row labels used throughout
# this project's scenario grids, see scenario_theta0 in
# code/sim_exp/bin/prior_comparison_sample.R and analogues) arranges its 3
# historical datasets' discrepancy from the current study.
#
# This is purely illustrative -- positions are HAND-CHOSEN to convey the
# qualitative pattern each Order label names, not computed from any
# simulated scenario's actual theta0/eps values (hence "(schematic)" on
# the x-axis, matching the reference figure this reproduces):
#   - Compatible:   discrepancy strictly decreases 1 -> 2 -> 3 (study 3 is
#                   closest to the current study, study 1 farthest) --
#                   exactly the order ONPP's eta1 <= eta2 <= eta3
#                   restriction assumes.
#   - Semi-compat.: close but not a strict staircase (2 closest, then 3,
#                   then 1) -- partially, not fully, order-compatible.
#   - Incompatible: discrepancy strictly INCREASES 1 -> 2 -> 3 (study 1
#                   closest, study 3 farthest) -- the order ONPP's
#                   restriction gets backwards, its worst case.
#   - Neutral:      all 3 studies sit at roughly the same (nonzero)
#                   discrepancy -- no order signal either way.
#   - Cancellation: two studies land on opposite sides of the current
#                   study by about the same amount, so a naive avg would
#                   "cancel out" to look deceptively close -- see
#                   code/sim_exp/cancellation_density_plots.R for the
#                   actual scenario this caricatures.
# =============================================================================

library(dplyr)
library(tidyverse)

eta_colors <- c("1" = "#7FB3E8", "2" = "#3B6FB6", "3" = "#0D2B52")

current_data_x <- 0.34

schematic_df <- tribble(
  ~category,       ~study, ~x,
  "Compatible",     1,      0.80,
  "Compatible",     2,      0.56,
  "Compatible",     3,      0.34,
  "Semi-compat.",   1,      0.80,
  "Semi-compat.",   3,      0.56,
  "Semi-compat.",   2,      0.34,
  "Incompatible",   3,      0.80,
  "Incompatible",   2,      0.56,
  "Incompatible",   1,      0.34,
  "Neutral",        3,      0.73,
  "Neutral",        2,      0.72,
  "Neutral",        1,      0.71,
  "Cancellation",   2,      0.48,
  "Cancellation",   3,      0.34,
  "Cancellation",   1,      0.2
)

# First factor level sits at the BOTTOM of a discrete ggplot y-axis, so
# list levels bottom-to-top to get the reference figure's top-to-bottom
# row order (Compatible first/top, Cancellation last/bottom).
category_levels <- c("Cancellation", "Neutral", "Incompatible", "Semi-compat.", "Compatible")
schematic_df$category <- factor(schematic_df$category, levels = category_levels)
schematic_df$study     <- factor(schematic_df$study, levels = c(1, 2, 3))

row_lines_df <- data.frame(category = factor(category_levels, levels = category_levels),
                            x_start = 0.06, x_end = 0.97)

p_order_schematic <- ggplot() +
  geom_segment(data = row_lines_df,
               aes(x = x_start, xend = x_end, y = category, yend = category),
               color = "grey75", linewidth = 0.6) +
  geom_vline(xintercept = current_data_x, color = "#D62728", linetype = "dashed", linewidth = 0.7) +
  annotate("text", x = current_data_x, y = length(category_levels) + 0.65,
           label = "current data", color = "#D62728", fontface = "bold", size = 5, hjust = 0.5) +
  geom_point(data = schematic_df, aes(x = x, y = category, fill = study),
             shape = 21, size = 11, color = "white", stroke = 0.6) +
  geom_text(data = schematic_df, aes(x = x, y = category, label = study),
            color = "white", fontface = "bold", size = 6) +
  scale_fill_manual(values = eta_colors, guide = "none") +
  annotate("segment", x = 0.04, xend = 0.99, y = 0.25, yend = 0.25,
           arrow = arrow(length = unit(0.25, "cm"), type = "closed"),
           linewidth = 0.6, color = "grey30") +
  annotate("text", x = 0.515, y = -0.15, label = "discrepancy of study k (schematic)",
           size = 5, color = "grey20") +
  coord_cartesian(xlim = c(0, 1), ylim = c(-0.3, length(category_levels) + 0.85), clip = "off") +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 15) +
  theme(
    panel.grid       = element_blank(),
    axis.text.x      = element_blank(),
    axis.ticks       = element_blank(),
    axis.text.y      = element_text(size = 16, color = "grey20", hjust = 1),
    plot.margin      = margin(t = 10, r = 15, b = 25, l = 10)
  )

print(p_order_schematic)
ggsave("results/figures/poster/order_discrepancy_schematic.pdf",
       p_order_schematic, width = 7, height = 4.5)
