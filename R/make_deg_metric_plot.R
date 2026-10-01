#' Create a differential expression metric plot
#'
#' Generate a bubble plot visualising differential expression metrics
#' across cell types or clusters.
#'
#' Each point represents a differentially expressed gene or feature.
#' The x-axis displays log fold-change values, while point size
#' represents statistical significance using
#' `-log10(adjusted p-value)`. Points are coloured according to
#' cell type or cluster identity.
#'
#' @param df A data frame or tibble containing differential expression
#'   results.
#' @param ct_col Character scalar specifying the column containing
#'   cell type, cluster, or spot annotations.
#' @param lfc_col Character scalar specifying the column containing
#'   log fold-change values.
#' @param padj_col Character scalar specifying the column containing
#'   adjusted p-values.
#' @param color_vec Optional named character vector of colours
#'   corresponding to cell types. If `NULL`, a default discrete
#'   palette is generated.
#' @param color_lab Optional character scalar specifying the colour
#'   legend title. Defaults to `ct_col`.
#' @param size_lab Optional character scalar specifying the point-size
#'   legend title. Defaults to `"-log10(p-adj)"`.
#' @param point_alpha Numeric scalar controlling point transparency.
#'   Default is `0.5`.
#' @param x_lab Character scalar specifying the x-axis label.
#'   Default is `"Log2 Fold Change"`.
#' @param y_lab Optional character scalar specifying the y-axis label.
#'   Default is `NULL`.
#' @param plot_title Optional character scalar specifying a plot title.
#' @param base_font_size Numeric scalar controlling the base font size.
#'   Default is `14`.
#'
#' @details
#' Cell type ordering follows factor levels of `ct_col` when supplied
#' as a factor. Otherwise, cell types are ordered alphabetically.
#'
#' Adjusted p-values less than or equal to zero are capped at
#' `1e-300` prior to transformation to avoid infinite values when
#' calculating `-log10(adjusted p-value)`.
#'
#' @return
#' A `ggplot2` object.
#'
#' @examples
#' make_deg_metric_plot(
#'   df = deg_results_df,
#'   ct_col = "cell_type",
#'   lfc_col = "log2FoldChange",
#'   padj_col = "padj",
#'   plot_title = "DEA Metrics Across Cell Types"
#' )
#'
#' @importFrom rlang .data
#' @importFrom ggtext element_textbox_simple
#'
#' @export
make_deg_metric_plot <- function(
    df,
    ct_col,
    lfc_col,
    padj_col,
    color_vec = NULL,
    color_lab = NULL,
    size_lab = NULL,
    point_alpha = 0.5,
    x_lab = "Log2 Fold Change",
    y_lab = NULL,
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
    c(
      ct_col,
      lfc_col,
      padj_col
    )
  )

  if (!is.numeric(df[[lfc_col]])) {
    stop(
      "Column '", lfc_col,
      "' must be numeric.",
      call. = FALSE
    )
  }

  if (!is.numeric(df[[padj_col]])) {
    stop(
      "Column '", padj_col,
      "' must be numeric.",
      call. = FALSE
    )
  }

  if (any(df[[padj_col]] < 0, na.rm = TRUE)) {
    stop(
      "Adjusted p-values must be non-negative.",
      call. = FALSE
    )
  }

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
  # Legend labels
  # -------------------------

  color_title <- if (
    is.null(color_lab)
  ) {
    ct_col
  } else {
    color_lab
  }

  size_title <- if (
    is.null(size_lab)
  ) {
    "-log10(p-adj)"
  } else {
    size_lab
  }

  # -------------------------
  # Process data
  # -------------------------

  df_processed <- df

  df_processed[[ct_col]] <- factor(
    df_processed[[ct_col]],
    levels = rev(ct_levels)
  )

  padj_values <- df_processed[[padj_col]]

  padj_values[
    padj_values <= 0 &
      !is.na(padj_values)
  ] <- 1e-300

  df_processed$minus_log10_padj <-
    -log10(padj_values)

  # -------------------------
  # Plot
  # -------------------------

  p <- ggplot2::ggplot(
    df_processed,
    ggplot2::aes(
      x = .data[[lfc_col]],
      y = .data[[ct_col]]
    )
  ) +
    ggplot2::geom_point(
      ggplot2::aes(
        colour = .data[[ct_col]],
        size = .data[["minus_log10_padj"]]
      ),
      alpha = point_alpha
    ) +
    ggplot2::scale_color_manual(
      values = color_vec,
      drop = FALSE
    ) +
    ggplot2::labs(
      x = x_lab,
      y = y_lab,
      colour = color_title,
      size = size_title,
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
      legend.title = ggplot2::element_text(
        size = base_font_size
      )
    ) +
    ggplot2::guides(
      colour = ggplot2::guide_legend(
        override.aes = list(
          alpha = 1,
          size = 4
        )
      )
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
