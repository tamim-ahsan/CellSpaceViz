#' Create a differential expression analysis log fold-change heat map
#'
#' Generate a heat map displaying log2 fold-change values across
#' conditions, cell types, or clusters for a selected set of genes.
#'
#' Significance annotations are displayed within heatmap tiles using
#' asterisks (`*`, `**`, `***`) according to user-specified adjusted
#' p-value thresholds.
#'
#' @param df A data frame or tibble containing differential expression
#'   analysis results.
#' @param lfc_col Character scalar specifying the column containing
#'   log2 fold-change values.
#' @param padj_col Character scalar specifying the column containing
#'   adjusted p-values.
#' @param gene_name_col Character scalar specifying the column
#'   containing gene identifiers or symbols.
#' @param genes_to_show Character vector specifying genes to include
#'   in the heat map.
#' @param condition_col Character scalar specifying the column
#'   containing conditions, cell types, or clusters.
#' @param sig_levels Numeric vector of adjusted p-value significance
#'   thresholds ordered from most stringent to least stringent.
#' @param scale_fill_option Character scalar specifying the Viridis
#'   colour palette.
#' @param scale_fill_name Character scalar specifying the fill
#'   legend title.
#' @param scale_fill_trans Character scalar specifying the fill-scale
#'   transformation.
#' @param x_lab Character scalar specifying the x-axis label.
#' @param y_lab Character scalar specifying the y-axis label.
#' @param plot_title Optional character scalar specifying a plot title.
#' @param base_font_size Numeric scalar controlling the base font size.
#'
#' @details
#' The order of genes follows the order supplied in
#' `genes_to_show`.
#'
#' Conditions follow factor levels of `condition_col` when supplied
#' as a factor. Otherwise, conditions are ordered alphabetically.
#'
#' Significance annotations are generated dynamically according to
#' the number of supplied thresholds. For example, using
#' `c(0.01, 0.05, 0.1)` results in:
#'
#' * `***` for adjusted p-values < 0.01
#' * `**` for adjusted p-values < 0.05
#' * `*` for adjusted p-values < 0.1
#'
#' @return
#' A `ggplot2` object.
#'
#' @examples
#' make_dea_lfc_heatmap(
#'   df = dea_heatmap_df,
#'   lfc_col = "log2FoldChange",
#'   padj_col = "padj",
#'   gene_name_col = "gene",
#'   genes_to_show = c(
#'     "gene1",
#'     "gene2",
#'     "gene3",
#'     "gene4"
#'   ),
#'   condition_col = "cell_type",
#'   plot_title = "Example DEA LFC Heatmap"
#' )
#'
#' @importFrom rlang .data
#' @importFrom ggtext element_textbox_simple
#'
#' @export
make_dea_lfc_heatmap <- function(
    df,
    lfc_col,
    padj_col,
    gene_name_col,
    genes_to_show,
    condition_col,
    sig_levels = c(
      0.01,
      0.05,
      0.1
    ),
    scale_fill_option = "magma",
    scale_fill_name = "Log2 Fold Change",
    scale_fill_trans = "identity",
    x_lab = "Condition",
    y_lab = "Gene",
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
      lfc_col,
      padj_col,
      gene_name_col,
      condition_col
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

  if (
    is.null(genes_to_show) ||
    length(genes_to_show) == 0
  ) {
    stop(
      "'genes_to_show' must contain at least one gene.",
      call. = FALSE
    )
  }

  if (
    !is.numeric(sig_levels) ||
    length(sig_levels) == 0
  ) {
    stop(
      "'sig_levels' must be a numeric vector.",
      call. = FALSE
    )
  }

  valid_viridis <- c(
    "magma",
    "inferno",
    "plasma",
    "viridis",
    "cividis",
    "rocket",
    "mako",
    "turbo"
  )

  if (!scale_fill_option %in% valid_viridis) {
    stop(
      "'scale_fill_option' must be one of: ",
      paste(
        valid_viridis,
        collapse = ", "
      ),
      call. = FALSE
    )
  }

  sig_levels <- sort(
    sig_levels,
    decreasing = FALSE
  )

  # -------------------------
  # Filter data
  # -------------------------

  df_filtered <- df |>
    dplyr::filter(
      .data[[gene_name_col]] %in%
        genes_to_show
    )

  if (nrow(df_filtered) == 0) {
    stop(
      "None of the genes in 'genes_to_show' were found in 'df'.",
      call. = FALSE
    )
  }

  # -------------------------
  # Significance labels
  # -------------------------

  n_sig <- length(sig_levels)

  padj_values <- df_filtered[[padj_col]]

  significance <- rep(
    NA_character_,
    length(padj_values)
  )

  for (i in seq_len(n_sig)) {

    stars <- paste(
      rep(
        "*",
        n_sig - i + 1
      ),
      collapse = ""
    )

    significance[
      padj_values < sig_levels[i] &
        is.na(significance)
    ] <- stars

  }

  df_filtered$significance <- significance

  # -------------------------
  # Factor ordering
  # -------------------------

  gene_levels <- rev(
    unique(
      genes_to_show
    )
  )

  df_filtered[[gene_name_col]] <- factor(
    df_filtered[[gene_name_col]],
    levels = gene_levels
  )

  cond_levels <- .get_factor_levels(
    df[[condition_col]]
  )

  df_filtered[[condition_col]] <- factor(
    df_filtered[[condition_col]],
    levels = cond_levels
  )

  lfc_midpoint <- mean(
    range(
      df_filtered[[lfc_col]],
      na.rm = TRUE
    )
  )

  # -------------------------
  # Plot
  # -------------------------

  p <- ggplot2::ggplot(
    df_filtered,
    ggplot2::aes(
      x = .data[[condition_col]],
      y = .data[[gene_name_col]],
      fill = .data[[lfc_col]]
    )
  ) +
    ggplot2::geom_tile(
      colour = "white",
      linewidth = 0.5
    ) +
    ggplot2::geom_text(
      ggplot2::aes(
        label = .data[["significance"]],
        colour = dplyr::if_else(
          .data[[lfc_col]] >
            lfc_midpoint,
          "black",
          "white"
        )
      ),
      size = base_font_size * 0.45,
      na.rm = TRUE,
      show.legend = FALSE
    ) +
    ggplot2::scale_color_identity() +
    ggplot2::scale_fill_viridis_c(
      option = scale_fill_option,
      name = scale_fill_name,
      transform = scale_fill_trans
    ) +
    ggplot2::labs(
      x = x_lab,
      y = y_lab,
      title = plot_title
    ) +
    ggplot2::theme_classic(
      base_size = base_font_size
    ) +
    ggplot2::theme(
      axis.line = ggplot2::element_blank(),
      panel.grid.major =
        ggplot2::element_blank(),
      axis.text.x =
        ggplot2::element_text(
          angle = 45,
          hjust = 1,
          vjust = 1,
          size = base_font_size
        ),
      axis.text.y =
        ggplot2::element_text(
          face = "italic",
          size = base_font_size
        ),
      axis.title =
        ggplot2::element_text(
          size = base_font_size + 2
        ),
      legend.title =
        ggplot2::element_text(
          size = base_font_size
        ),
      legend.position = "right"
    )

  if (!is.null(plot_title)) {

    p <- p +
      ggplot2::theme(
        plot.title =
          ggtext::element_textbox_simple(
            size = base_font_size + 2,
            face = "bold",
            hjust = 0.5,
            halign = 0.5,
            linetype = 1,
            box.color = "black",
            linewidth = 0.8,
            padding =
              ggplot2::margin(
                5, 5, 5, 5
              ),
            margin =
              ggplot2::margin(
                b = 10
              )
          )
      )

  }

  return(p)

}
