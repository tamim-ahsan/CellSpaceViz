#' Create a GSEA pathway enrichment plot
#'
#' Generate a horizontal bar plot showing the most significantly
#' enriched pathways from Gene Set Enrichment Analysis (GSEA).
#'
#' Pathways are ordered according to normalized enrichment score (NES)
#' and coloured according to statistical significance using
#' `-log10(adjusted p-value)`.
#'
#' Positive enrichment scores represent positively enriched
#' pathways, whereas negative enrichment scores represent
#' negatively enriched pathways.
#'
#' The top positively enriched and negatively enriched pathways
#' are displayed separately according to the specified selection
#' criteria.
#'
#' @param df A data frame or tibble containing GSEA enrichment
#'   results.
#' @param enrichment_score_col Character scalar specifying the
#'   column containing enrichment scores (e.g. NES).
#' @param padj_col Character scalar specifying the column
#'   containing adjusted p-values.
#' @param description_col Character scalar specifying the column
#'   containing pathway names or descriptions.
#' @param ontology_col Optional character scalar specifying the
#'   ontology/category column.
#' @param ontology Optional character scalar specifying an ontology
#'   to retain when `ontology_col` is supplied.
#' @param n_pathways_to_show Integer scalar specifying the maximum
#'   number of positively enriched and negatively enriched pathways
#'   to display.
#' @param padj_threshold Numeric scalar specifying the adjusted
#'   p-value filtering threshold.
#' @param fill_scale_option Character scalar specifying the Viridis
#'   colour palette.
#' @param fill_scale_name Character scalar specifying the fill-scale
#'   legend title.
#' @param x_lab Character scalar specifying the x-axis label.
#' @param y_lab Optional character scalar specifying the y-axis label.
#' @param plot_title Optional character scalar specifying a plot title.
#' @param description_text_size Numeric scalar controlling pathway
#'   label size.
#' @param base_font_size Numeric scalar controlling the base font size.
#'
#' @details
#' Pathways are filtered according to `padj_threshold`.
#'
#' Positive pathways are ranked from the largest positive
#' enrichment score, whereas negative pathways are ranked from the
#' smallest negative enrichment score.
#'
#' Colours correspond to:
#'
#' `-log10(adjusted p-value)`
#'
#' so more statistically significant pathways appear with stronger
#' colour intensity.
#'
#' @return
#' A `ggplot2` object.
#'
#' @examples
#' make_gsea_pathway_plot(
#'   df = gsea_results_df,
#'   enrichment_score_col = "NES",
#'   padj_col = "padj",
#'   description_col = "Description",
#'   n_pathways_to_show = 5,
#'   plot_title = "Example GSEA Pathways"
#' )
#'
#' @importFrom rlang .data
#' @importFrom ggtext element_textbox_simple
#'
#' @export
#'

make_gsea_pathway_plot <- function(
    df,
    enrichment_score_col,
    padj_col,
    description_col,
    ontology_col = NULL,
    ontology = NULL,
    n_pathways_to_show = 10,
    padj_threshold = 0.05,
    fill_scale_option = "magma",
    fill_scale_name = "-log10(p-adj)",
    x_lab = "Normalized Enrichment Score (NES)",
    y_lab = NULL,
    plot_title = NULL,
    description_text_size = 4,
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
      enrichment_score_col,
      padj_col,
      description_col
    )
  )

  if (!is.numeric(df[[enrichment_score_col]])) {
    stop(
      "Column '", enrichment_score_col,
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

  if (
    !is.numeric(padj_threshold) ||
    length(padj_threshold) != 1 ||
    padj_threshold <= 0 ||
    padj_threshold >= 1
  ) {
    stop(
      "'padj_threshold' must be a single value between 0 and 1.",
      call. = FALSE
    )
  }

  if (
    !is.numeric(n_pathways_to_show) ||
    n_pathways_to_show < 1
  ) {
    stop(
      "'n_pathways_to_show' must be a positive integer.",
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

  if (!fill_scale_option %in% valid_viridis) {
    stop(
      "'fill_scale_option' must be one of: ",
      paste(valid_viridis, collapse = ", "),
      call. = FALSE
    )
  }

  # -------------------------
  # Ontology filtering
  # -------------------------

  if (!is.null(ontology_col) ||
      !is.null(ontology)) {

    if (
      is.null(ontology_col) ||
      is.null(ontology)
    ) {

      stop(
        "Both 'ontology_col' and 'ontology' must be supplied.",
        call. = FALSE
      )

    }

    .validate_columns(
      df,
      ontology_col
    )

    df <- dplyr::filter(
      df,
      .data[[ontology_col]] ==
        ontology
    )

  }

  # -------------------------
  # Significant pathways
  # -------------------------

  df_sig <- dplyr::filter(
    df,
    .data[[padj_col]] <
      padj_threshold
  )

  if (nrow(df_sig) == 0) {
    stop(
      "No pathways passed the adjusted p-value threshold.",
      call. = FALSE
    )
  }

  df_up <- df_sig |>
    dplyr::filter(
      .data[[enrichment_score_col]] > 0
    ) |>
    dplyr::arrange(
      dplyr::desc(
        .data[[enrichment_score_col]]
      )
    ) |>
    dplyr::slice_head(
      n = n_pathways_to_show
    )

  df_down <- df_sig |>
    dplyr::filter(
      .data[[enrichment_score_col]] < 0
    ) |>
    dplyr::arrange(
      .data[[enrichment_score_col]]
    ) |>
    dplyr::slice_head(
      n = n_pathways_to_show
    )

  df_final <- dplyr::bind_rows(
    df_up,
    df_down
  )

  if (nrow(df_final) == 0) {
    stop(
      "No pathways available after filtering.",
      call. = FALSE
    )
  }

  padj_values <- df_final[[padj_col]]

  padj_values[
    padj_values <= 0 &
      !is.na(padj_values)
  ] <- 1e-300

  df_final$minus_log10_padj <-
    -log10(padj_values)

  max_nes <- max(
    df_final[[enrichment_score_col]],
    na.rm = TRUE
  )

  min_nes <- min(
    df_final[[enrichment_score_col]],
    na.rm = TRUE
  )

  x_limits <- c(
    ifelse(
      min_nes < 0,
      min_nes * 1.2,
      -0.1 * max_nes
    ),
    ifelse(
      max_nes > 0,
      max_nes * 1.2,
      -0.1 * min_nes
    )
  )

  # -------------------------
  # Plot
  # -------------------------

  p <- ggplot2::ggplot(
    df_final,
    ggplot2::aes(
      x = .data[[enrichment_score_col]],
      y = stats::reorder(
        .data[[description_col]],
        .data[[enrichment_score_col]]
      ),
      fill = .data[["minus_log10_padj"]]
    )
  ) +
    ggplot2::geom_col(
      width = 0.8,
      colour = "black",
      linewidth = 0.5
    ) +
    ggplot2::scale_fill_viridis_c(
      option = fill_scale_option,
      name = fill_scale_name
    ) +
    ggplot2::coord_cartesian(
      xlim = x_limits,
      clip = "off"
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
      axis.title = ggplot2::element_text(
        size = base_font_size + 2
      ),
      axis.text.x = ggplot2::element_text(
        size = base_font_size
      ),
      legend.title = ggplot2::element_text(
        size = base_font_size
      ),
      plot.margin = ggplot2::margin(
        1,
        1,
        1,
        1,
        "cm"
      )
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
