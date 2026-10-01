#' Create a cell, nucleus, or spot count bar plot
#'
#' Generate a horizontal bar plot showing the total number of cells,
#' nuclei, or spots per cell type or cluster.
#'
#' Counts are automatically calculated for each unique group in
#' `ct_col`. Text labels displaying exact counts are rendered at the
#' end of each bar.
#'
#' @param df A data frame or tibble containing cell metadata.
#' @param ct_col Character scalar specifying the column containing
#'   cell type, cluster, or spot annotations.
#' @param color_vec Optional named character vector of colours
#'   corresponding to cell types. If `NULL`, a default discrete
#'   palette is generated.
#' @param label_size Numeric scalar controlling the size of count
#'   labels. Default is `3.5`.
#' @param cell_or_spot Character scalar specifying whether the
#'   observations represent `"Cell"`, `"Nucleus"`, or `"Spot"`.
#' @param plot_title Optional character scalar specifying a plot
#'   title.
#' @param base_font_size Numeric scalar controlling the base font
#'   size used throughout the plot.
#'
#' @details
#' Cell type ordering follows the factor levels of `ct_col` when
#' supplied as a factor. Otherwise, cell types are ordered
#' alphabetically. The first factor level is displayed at the top
#' of the figure.
#'
#' @return
#' A `ggplot2` object.
#'
#' @examples
#' make_nuclei_count_barplot(
#'   df = nuclei_count_df,
#'   ct_col = "cell_type",
#'   cell_or_spot = "Nucleus",
#'   plot_title = "Example Nuclei Counts"
#' )
#' @importFrom rlang .data
#'
#' @export
make_nuclei_count_barplot <- function(
    df,
    ct_col,
    color_vec = NULL,
    label_size = 3.5,
    cell_or_spot = "Cell",
    plot_title = NULL,
    base_font_size = 14
) {

  # -------------------------
  # Input validation
  # -------------------------

  if (!inherits(df, "data.frame")) {
    stop(
      "'df' must be a data.frame or tibble.",
      call. = FALSE
    )
  }

  .validate_columns(
    df,
    ct_col
  )

  cell_or_spot <- match.arg(
    cell_or_spot,
    choices = c(
      "Cell",
      "Nucleus",
      "Spot"
    )
  )

  # -------------------------
  # Cell-type handling
  # -------------------------

  ct_levels <- .get_factor_levels(
    df[[ct_col]]
  )

  if (length(ct_levels) == 0) {
    stop(
      "No valid (non-NA) values found in 'ct_col'.",
      call. = FALSE
    )
  }

  color_vec <- .validate_colors(
    levels = ct_levels,
    color_vec = color_vec
  )

  # -------------------------
  # Process data
  # -------------------------

  counts_table <- table(
    factor(
      df[[ct_col]],
      levels = ct_levels
    )
  )

  plot_df <- data.frame(
    celltype = factor(
      names(counts_table),
      levels = rev(ct_levels)
    ),
    count = as.numeric(counts_table),
    stringsAsFactors = FALSE
  )

  plot_df$label_text <- format(
    plot_df$count,
    big.mark = ",",
    scientific = FALSE
  )

  # -------------------------
  # Plot
  # -------------------------

  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(
      x = .data[["count"]],
      y = .data[["celltype"]],
      fill = .data[["celltype"]],
      label = .data[["label_text"]]
    )
  ) +
    ggplot2::geom_col(
      show.legend = FALSE
    ) +
    ggplot2::geom_text(
      hjust = -0.15,
      size = label_size
    ) +
    ggplot2::scale_fill_manual(
      values = color_vec,
      drop = FALSE
    ) +
    ggplot2::scale_x_continuous(
      expand = ggplot2::expansion(
        mult = c(0, 0.15)
      )
    ) +
    ggplot2::coord_cartesian(
      clip = "off"
    ) +
    ggplot2::labs(
      x = paste(cell_or_spot, "Count"),
      y = NULL,
      title = plot_title
    ) +
    ggplot2::theme_classic(
      base_size = base_font_size
    ) +
    ggplot2::theme(
      axis.title = ggplot2::element_text(
        size = base_font_size + 2
      ),
      axis.text.y = ggplot2::element_text(
        size = base_font_size
      ),
      legend.position = "none"
    )

  if (!is.null(plot_title)) {

    p <- p +
      ggplot2::theme(
        plot.title = ggtext::element_textbox_simple(
          size = base_font_size + 2,
          face = "bold",
          hjust = 0.5,
          halign = 0.5,
          linetype = 1,
          box.color = "black",
          linewidth = 0.8,
          width = NULL,
          padding = ggplot2::margin(
            5, 5, 5, 5
          ),
          margin = ggplot2::margin(
            b = 10
          )
        )
      )

  }

  return(p)

}
