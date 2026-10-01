#' Create a differential expression analysis volcano plot
#'
#' Generate a volcano plot displaying differential expression analysis
#' results using log2 fold-change and adjusted p-value statistics.
#'
#' Each point represents a gene or feature. The x-axis displays
#' log2 fold-change values and the y-axis displays
#' `-log10(adjusted p-value)`.
#'
#' Features are classified as up-regulated, down-regulated, or
#' insignificant according to the specified significance threshold.
#' The most significant up-regulated and down-regulated genes can
#' optionally be labelled using `ggrepel`.
#'
#' @param df A data frame or tibble containing differential expression
#'   results.
#' @param lfc_col Character scalar specifying the column containing
#'   log2 fold-change values.
#' @param padj_col Character scalar specifying the column containing
#'   adjusted p-values.
#' @param gene_label_col Character scalar specifying the column
#'   containing gene identifiers or symbols.
#' @param color_vec Optional named character vector specifying colours
#'   for `"Upregulated"`, `"Downregulated"`, and `"Insignificant"`
#'   categories.
#' @param x_lab Character scalar specifying the x-axis label.
#' @param y_lab Character scalar specifying the y-axis label.
#' @param color_lab Character scalar specifying the colour legend title.
#' @param plot_title Optional character scalar specifying a plot title.
#' @param padj_thresh Numeric scalar specifying the adjusted p-value
#'   significance threshold.
#' @param point_alpha Numeric scalar controlling point transparency.
#' @param point_size Numeric scalar controlling point size.
#' @param base_font_size Numeric scalar controlling the base font size.
#' @param gene_text_size Numeric scalar controlling gene-label size.
#' @param sig_gene_to_show Integer scalar specifying the number of
#'   up-regulated and down-regulated genes to label.
#'
#' @details
#' Adjusted p-values less than or equal to zero are capped at
#' `1e-300` before transformation to avoid infinite values when
#' calculating `-log10(adjusted p-value)`.
#'
#' Up-regulated genes satisfy:
#'
#' * `padj < padj_thresh`
#' * `log2 fold-change > 0`
#'
#' Down-regulated genes satisfy:
#'
#' * `padj < padj_thresh`
#' * `log2 fold-change < 0`
#'
#' All remaining genes are classified as insignificant.
#'
#' @return
#' A `ggplot2` object.
#'
#' @examples
#' make_dea_volcano_plot(
#'   df = dea_volcano_df,
#'   lfc_col = "log2FoldChange",
#'   padj_col = "padj",
#'   gene_label_col = "gene",
#'   padj_thresh = 0.05,
#'   plot_title = "Example DEA Volcano Plot"
#' )
#'
#' @importFrom rlang .data
#' @importFrom ggrepel geom_label_repel
#' @importFrom ggtext element_textbox_simple
#'
#' @export
make_dea_volcano_plot <- function(
    df,
    lfc_col,
    padj_col,
    gene_label_col,
    color_vec = NULL,
    x_lab = "Log2 Fold Change",
    y_lab = "-log10(p-adj)",
    color_lab = "Significance",
    plot_title = NULL,
    padj_thresh = 0.1,
    point_alpha = 0.5,
    point_size = 4,
    base_font_size = 14,
    gene_text_size = 4,
    sig_gene_to_show = 5
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
      gene_label_col
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

  if (
    !is.numeric(padj_thresh) ||
    length(padj_thresh) != 1 ||
    padj_thresh <= 0 ||
    padj_thresh >= 1
  ) {
    stop(
      "'padj_thresh' must be a single value between 0 and 1.",
      call. = FALSE
    )
  }

  if (
    !is.numeric(sig_gene_to_show) ||
    length(sig_gene_to_show) != 1 ||
    sig_gene_to_show < 0
  ) {
    stop(
      "'sig_gene_to_show' must be a non-negative integer.",
      call. = FALSE
    )
  }

  # -------------------------
  # Colours
  # -------------------------

  categories <- c(
    "Upregulated",
    "Downregulated",
    "Insignificant"
  )

  if (is.null(color_vec)) {

    color_vec <- c(
      "Upregulated" = "#E41A1C",
      "Downregulated" = "#377EB8",
      "Insignificant" = "grey70"
    )

  } else {

    if (length(color_vec) != 3) {
      stop(
        "'color_vec' must contain exactly 3 colours.",
        call. = FALSE
      )
    }

    if (is.null(names(color_vec))) {

      names(color_vec) <- categories

    } else {

      missing_categories <- setdiff(
        categories,
        names(color_vec)
      )

      if (length(missing_categories) > 0) {

        stop(
          "Missing colour mapping for category/categories: ",
          paste(
            missing_categories,
            collapse = ", "
          ),
          call. = FALSE
        )

      }

      color_vec <- color_vec[categories]

    }

  }

  # -------------------------
  # Process data
  # -------------------------

  padj_values <- df[[padj_col]]

  padj_values[
    padj_values <= 0 &
      !is.na(padj_values)
  ] <- 1e-300

  df_processed <- df |>
    dplyr::mutate(
      minus_log10_padj = -log10(padj_values),
      significance = dplyr::case_when(
        is.na(.data[[padj_col]]) ~ "Insignificant",
        .data[[padj_col]] < padj_thresh &
          .data[[lfc_col]] > 0 ~ "Upregulated",
        .data[[padj_col]] < padj_thresh &
          .data[[lfc_col]] < 0 ~ "Downregulated",
        TRUE ~ "Insignificant"
      ),
      significance = factor(
        .data[["significance"]],
        levels = categories
      )
    )

  # -------------------------
  # Gene labels
  # -------------------------

  up_count <- sum(
    df_processed$significance == "Upregulated",
    na.rm = TRUE
  )

  down_count <- sum(
    df_processed$significance == "Downregulated",
    na.rm = TRUE
  )

  if (
    up_count < sig_gene_to_show ||
    down_count < sig_gene_to_show
  ) {

    warning(
      "Number of significant up-regulated (",
      up_count,
      ") or down-regulated (",
      down_count,
      ") genes is less than 'sig_gene_to_show' (",
      sig_gene_to_show,
      ").",
      call. = FALSE
    )

  }

  df_processed <- df_processed |>
    dplyr::arrange(
      dplyr::desc(
        .data[["minus_log10_padj"]]
      )
    ) |>
    dplyr::mutate(
      is_up = .data[["significance"]] ==
        "Upregulated",
      is_down = .data[["significance"]] ==
        "Downregulated",
      rank_up = cumsum(
        .data[["is_up"]]
      ),
      rank_down = cumsum(
        .data[["is_down"]]
      ),
      label_gene = dplyr::case_when(
        .data[["is_up"]] &
          .data[["rank_up"]] <= sig_gene_to_show ~ TRUE,
        .data[["is_down"]] &
          .data[["rank_down"]] <= sig_gene_to_show ~ TRUE,
        TRUE ~ FALSE
      )
    )

  # -------------------------
  # Plot
  # -------------------------

  p <- ggplot2::ggplot(
    df_processed,
    ggplot2::aes(
      x = .data[[lfc_col]],
      y = .data[["minus_log10_padj"]],
      colour = .data[["significance"]]
    )
  ) +
    ggplot2::geom_point(
      alpha = point_alpha,
      size = point_size
    ) +
    ggplot2::geom_hline(
      yintercept = -log10(padj_thresh),
      linetype = "dashed",
      colour = "grey50"
    ) +
    ggrepel::geom_label_repel(
      data = dplyr::filter(
        df_processed,
        .data[["label_gene"]]
      ),
      ggplot2::aes(
        label = .data[[gene_label_col]]
      ),
      fontface = "italic",
      max.overlaps = Inf,
      size = gene_text_size,
      show.legend = FALSE
    ) +
    ggplot2::scale_color_manual(
      values = color_vec,
      drop = FALSE
    ) +
    ggplot2::labs(
      x = x_lab,
      y = y_lab,
      colour = color_lab,
      title = plot_title
    ) +
    ggplot2::theme_classic(
      base_size = base_font_size
    ) +
    ggplot2::theme(
      axis.title = ggplot2::element_text(
        size = base_font_size + 2
      ),
      axis.text = ggplot2::element_text(
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
          size = point_size
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
