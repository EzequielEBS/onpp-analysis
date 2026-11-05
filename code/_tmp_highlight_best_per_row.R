
# =============================================================================
# Highlight only the best (smallest) value per row in a gt summary table
# =============================================================================
#
# Replaces the old per-row data_color() gradient (used in analysis_bin.r,
# analysis_normal.r, analysis_lm.R's build_summary_gt_table()) that colored
# every cell in a row on a continuous mako/plasma scale and let gt
# auto-pick black-or-white text per cell for contrast. With cells often
# landing at similar shades, it wasn't obvious at a glance which cell (if
# any) was actually the best method for that row, and the black/white text
# felt arbitrary. This instead marks only the row's minimum with a solid
# fill + bold white text and leaves every other cell plain (white
# background, black text) -- exactly one highlighted cell per row/group,
# everything else uniform.
#
# Arguments:
#   tbl     a gt table object, already built (gt(), cols_label(), spanners,
#           etc. all applied)
#   df      the data.frame the table was built from -- row order must match
#           `tbl` (i.e. no filtering/reordering happened in between)
#   columns character vector of column names that participate in the
#           highlighting (typically the method columns, e.g. c("NPP",
#           "NPP_SEQ", "ONPP", "ONPP_SEQ")) -- also the columns reset to a
#           plain white/black style before highlighting is applied, so any
#           earlier gt_color_rows()/data_color() coloring on them is wiped
#           out first.
#   groups  optional list of character vectors partitioning `columns` into
#           groups that should each get their own "best" cell per row,
#           instead of one best cell across all of `columns`. Use this
#           when `columns` mixes quantities on different scales that
#           shouldn't be compared against each other -- e.g. lm's summary
#           tables hold both beta and sigma columns per row, so the best
#           method should be picked separately within the beta columns and
#           within the sigma columns: groups = list(beta_cols, sigma_cols).
#           Defaults to list(columns), i.e. one group covering everything.
#   color       fill color for the highlighted (best) cell
#   text_color  text color for the highlighted cell (default white, for
#               contrast against the default green `color`)
#
# Returns the gt table with the highlighting applied.
#
# Example:
#   table_mse_theta <- highlight_best_per_row(table_mse_theta, df_mse_theta,
#                                              columns = method_cols,
#                                              color = "#2E7D32")  # MSE/norm: green
#   table_wis_theta <- highlight_best_per_row(table_wis_theta, df_wis_theta,
#                                              columns = method_cols,
#                                              color = "#6A1B9A")  # WIS: purple
#   tbl <- highlight_best_per_row(tbl, df, columns = method_cols,
#                                  groups = list(beta_cols, sigma_cols),
#                                  color = "#2E7D32")
highlight_best_per_row <- function(tbl, df, columns, groups = NULL,
                                    color = "#2E7D32", text_color = "white") {
  if (is.null(groups)) groups <- list(columns)

  # Reset every candidate cell to a plain, uniform style first so no
  # leftover coloring from an earlier gt_color_rows()/data_color() call
  # survives on the cells that don't end up highlighted.
  tbl <- tbl |>
    tab_style(
      style     = list(cell_fill(color = "white"), cell_text(color = "black")),
      locations = cells_body(columns = all_of(columns))
    )

  for (i in seq_len(nrow(df))) {
    for (grp in groups) {
      row_vals <- as.numeric(df[i, grp])
      best_col <- grp[which.min(row_vals)]
      tbl <- tbl |>
        tab_style(
          style = list(cell_fill(color = color),
                      cell_text(color = text_color, weight = "bold")),
          locations = cells_body(columns = all_of(best_col), rows = i)
        )
    }
  }
  tbl
}
