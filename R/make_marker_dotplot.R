#' Create a marker-gene dot plot
#'
#' Generate a dot plot showing average gene expression and the
#' percentage of cells or spots expressing selected marker genes
#' across cell types or clusters.
#'
#' Dot color represents marker gene expression and dot size represents
#' the percentage of cells or spots with non-zero expression.
#' Marker genes can be grouped into categories that are displayed
#' as separate facet panels.
#'
#' @param df A data frame or tibble containing cell metadata and
#'   expression values.
#' @param marker_genes Character vector specifying marker genes
#'   to display. Each gene must correspond to a numeric column in
#'   `df`.
#' @param ct_col Character scalar specifying the column containing
#'   cell type, cluster, or spot annotations.
#' @param gene_groups Character vector assigning each marker gene
#'   to a gene group for faceted visualization. Must have the same
#'   length as `marker_genes`.
#' @param color_scale_option Character scalar specifying the
#'   viridis color palette. One of `"magma"`, `"inferno"`,
#'   `"plasma"`, `"viridis"`, `"cividis"`, `"rocket"`,
#'   `"mako"`, or `"turbo"`.
#' @param color_scale_name Optional character scalar specifying
#'   the legend title for the expression color scale.
#' @param plot_title Optional character scalar specifying a plot
#'   title.
#' @param cells_or_spots Character scalar specifying whether the
#'   observations represent `"cells"` or `"spots"`.
#' @param base_font_size Numeric scalar controlling the base font
#'   size used throughout the plot.
#'
#' @details
#' For each marker gene and cell type combination, the function
#' calculates:
#'
#' * Mean expression (`mean_count`)
#' * Percentage of cells or spots with expression greater than zero
#'
#' Mean expression is displayed using color and the percentage of
#' expressing cells or spots is displayed using point size.
#'
#' Cell type ordering follows the factor levels of `ct_col` when
#' supplied as a factor. Otherwise, cell types are ordered
#' alphabetically.
#'
#' @return
#' A `ggplot2` object.
#'
#' @examples
#' marker_genes <- c(
#'   "geneA1", "geneA2", "geneA3",
#'   "geneB1", "geneB2", "geneB3",
#'   "geneC1", "geneC2", "geneC3",
#'   "geneD1", "geneD2", "geneD3"
#' )
#'
#' gene_groups <- c(
#'   rep("Cell-type A", 3),
#'   rep("Cell-type B", 3),
#'   rep("Cell-type C", 3),
#'   rep("Cell-type D", 3)
#' )
#'
#' make_marker_dotplot(
#'   df = marker_df,
#'   marker_genes = marker_genes,
#'   ct_col = "cell_type",
#'   gene_groups = gene_groups,
#'   plot_title = "Example Marker Dot Plot"
#' )
#'
#' @importFrom rlang .data :=
#'
#' @export

make_marker_dotplot <- function(
    df,
    marker_genes,
    ct_col,
    gene_groups,
    color_scale_option = "magma",
    color_scale_name = NULL,
    plot_title = NULL,
    cells_or_spots = "cells",
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

  if (!ct_col %in% colnames(df)) {
    stop(
      "Column '", ct_col,
      "' not found in 'df'.",
      call. = FALSE
    )
  }

  if (length(marker_genes) == 0) {
    stop(
      "'marker_genes' must contain at least one gene.",
      call. = FALSE
    )
  }

  missing_genes <- setdiff(
    marker_genes,
    colnames(df)
  )

  if (length(missing_genes) > 0) {
    stop(
      "Missing marker gene column(s): ",
      paste(missing_genes, collapse = ", "),
      call. = FALSE
    )
  }

  if (anyDuplicated(marker_genes)) {
    stop(
      "'marker_genes' contains duplicated genes.",
      call. = FALSE
    )
  }

  if (length(marker_genes) != length(gene_groups)) {
    stop(
      "Length of 'marker_genes' (",
      length(marker_genes),
      ") does not match length of 'gene_groups' (",
      length(gene_groups),
      ").",
      call. = FALSE
    )
  }

  if (anyNA(gene_groups)) {
    stop(
      "'gene_groups' contains missing values.",
      call. = FALSE
    )
  }

  non_numeric <- marker_genes[
    !vapply(
      df[marker_genes],
      is.numeric,
      logical(1)
    )
  ]

  if (length(non_numeric) > 0) {
    stop(
      "Marker gene columns must be numeric: ",
      paste(non_numeric, collapse = ", "),
      call. = FALSE
    )
  }

  cells_or_spots <- match.arg(
    cells_or_spots,
    choices = c(
      "cells",
      "spots"
    )
  )

  valid_scales <- c(
    "magma",
    "inferno",
    "plasma",
    "viridis",
    "cividis",
    "rocket",
    "mako",
    "turbo"
  )

  if (!color_scale_option %in% valid_scales) {
    stop(
      "'color_scale_option' must be one of: ",
      paste(valid_scales, collapse = ", "),
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

  ct_levels_y <- rev(ct_levels)

  # -------------------------
  # Gene grouping
  # -------------------------

  gene_group_map <- data.frame(
    genes = marker_genes,
    gene_groups = factor(
      gene_groups,
      levels = unique(gene_groups)
    ),
    stringsAsFactors = FALSE
  )

  # -------------------------
  # Process data
  # -------------------------

  df_processed <- df |>
    tidyr::pivot_longer(
      cols = tidyselect::all_of(marker_genes),
      names_to = "genes",
      values_to = "counts"
    )

  df_processed[[ct_col]] <- factor(
    df_processed[[ct_col]],
    levels = ct_levels_y
  )

  df_processed <- df_processed |>
    dplyr::mutate(
      genes = factor(
        .data[["genes"]],
        levels = marker_genes
      )
    ) |>
    dplyr::group_by(
      .data[[ct_col]],
      .data[["genes"]]
    ) |>
    dplyr::summarise(
      mean_expression = mean(
        .data[["counts"]],
        na.rm = TRUE
      ),
      pct_expressing = mean(
        .data[["counts"]] > 0,
        na.rm = TRUE
      ) * 100,
      .groups = "drop"
    ) |>
    dplyr::left_join(
      gene_group_map,
      by = "genes"
    )

  # -------------------------
  # Labels
  # -------------------------

  color_title <- if (is.null(color_scale_name)) {
    "Mean\nExpression"
  } else {
    color_scale_name
  }

  size_title <- paste0(
    "Percentage of ",
    cells_or_spots,
    "\nexpressing the gene"
  )

  # -------------------------
  # Plot
  # -------------------------

  p <- ggplot2::ggplot(
    df_processed,
    ggplot2::aes(
      x = .data[["genes"]],
      y = .data[[ct_col]],
      colour = .data[["mean_expression"]],
      size = .data[["pct_expressing"]]
    )
  ) +
    ggplot2::geom_point() +
    ggplot2::facet_grid(
      . ~ .data[["gene_groups"]],
      scales = "free_x",
      space = "free_x"
    ) +
    ggplot2::scale_color_viridis_c(
      option = color_scale_option,
      name = color_title
    ) +
    ggplot2::scale_size_continuous(
      name = size_title
    ) +
    ggplot2::labs(
      x = "Marker Genes",
      y = NULL,
      title = plot_title
    ) +
    ggplot2::theme_classic(
      base_size = base_font_size
    ) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        angle = 45,
        hjust = 1,
        vjust = 1
      ),
      axis.title = ggplot2::element_text(
        size = base_font_size + 2
      ),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(
        face = "bold",
        size = base_font_size + 1
      ),
      panel.spacing.x = ggplot2::unit(
        1.5,
        "lines"
      ),
      panel.border = ggplot2::element_rect(
        colour = "grey80",
        fill = NA
      ),
      legend.title = ggplot2::element_text(
        size = base_font_size
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
