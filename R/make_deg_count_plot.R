#' Create a differentially expressed gene (DEG) count bar plot
#'
#' Generate a stacked horizontal bar plot showing the numbers of
#' up-regulated and down-regulated differentially expressed genes
#' (DEGs) across cell types or clusters.
#'
#' Bars are divided into up-regulated and down-regulated gene counts.
#' Counts for each segment are displayed within the bar, and total DEG
#' counts are displayed at the end of each bar.
#'
#' @param df A data frame or tibble containing DEG summary statistics.
#' @param ct_col Character scalar specifying the column containing
#'   cell type or cluster annotations.
#' @param deg_num_col Character scalar specifying the column containing
#'   total DEG counts.
#' @param up_col Character scalar specifying the column containing
#'   counts of up-regulated genes.
#' @param down_col Character scalar specifying the column containing
#'   counts of down-regulated genes.
#' @param text_size Numeric scalar controlling the size of count labels.
#' @param stack_text_color Character scalar specifying the colour of
#'   labels displayed within stacked bar segments.
#' @param x_lab Character scalar specifying the x-axis label.
#' @param y_lab Optional character scalar specifying the y-axis label.
#' @param plot_title Optional character scalar specifying a plot title.
#' @param color_vec Optional named character vector containing colours
#'   for up-regulated and down-regulated categories.
#' @param base_font_size Numeric scalar controlling the base font size.
#'
#' @details
#' Cell type ordering follows factor levels of `ct_col` when supplied
#' as a factor. Otherwise, cell types are ordered alphabetically.
#'
#' @return
#' A `ggplot2` object.
#'
#' @examples
#' make_deg_count_plot(
#'   df = deg_summary_df,
#'   ct_col = "cell_type",
#'   deg_num_col = "total_degs",
#'   up_col = "up_regulated",
#'   down_col = "down_regulated"
#' )
#'
#' @importFrom rlang .data
#'
#' @export

make_deg_count_plot <- function(
    df,
    ct_col,
    deg_num_col,
    up_col,
    down_col,
    text_size = 4,
    stack_text_color = "white",
    x_lab = "Number of DEGs",
    y_lab = NULL,
    plot_title = NULL,
    color_vec = NULL,
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
      deg_num_col,
      up_col,
      down_col
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

  ct_levels_y <- rev(ct_levels)

  # -------------------------
  # Colours
  # -------------------------

  if (is.null(color_vec)) {

    color_vec <- c(
      "#E41A1C",
      "#377EB8"
    )

    names(color_vec) <- c(
      up_col,
      down_col
    )

  } else {

    if (length(color_vec) != 2) {
      stop(
        "'color_vec' must contain exactly 2 colours.",
        call. = FALSE
      )
    }

    if (is.null(names(color_vec))) {

      names(color_vec) <- c(
        up_col,
        down_col
      )

    } else {

      missing_groups <- setdiff(
        c(up_col, down_col),
        names(color_vec)
      )

      if (length(missing_groups) > 0) {
        stop(
          "Missing colour mapping for: ",
          paste(missing_groups, collapse = ", "),
          call. = FALSE
        )
      }

      color_vec <- color_vec[
        c(up_col, down_col)
      ]

    }

  }

  # -------------------------
  # Data processing
  # -------------------------

  total_label_df <- data.frame(
    celltype = factor(
      df[[ct_col]],
      levels = ct_levels_y
    ),
    total_deg = df[[deg_num_col]],
    label_text = format(
      df[[deg_num_col]],
      big.mark = ",",
      scientific = FALSE
    ),
    stringsAsFactors = FALSE
  )

  df_processed <- df |>
    tidyr::pivot_longer(
      cols = c(
        tidyselect::all_of(up_col),
        tidyselect::all_of(down_col)
      ),
      names_to = "regulation",
      values_to = "deg_count"
    ) |>
    dplyr::mutate(
      celltype = factor(
        .data[[ct_col]],
        levels = ct_levels_y
      ),
      regulation = factor(
        .data[["regulation"]],
        levels = c(up_col, down_col)
      ),
      segment_label = ifelse(
        .data[["deg_count"]] > 0,
        format(
          .data[["deg_count"]],
          big.mark = ",",
          scientific = FALSE
        ),
        ""
      )
    )

  # -------------------------
  # Plot
  # -------------------------

  p <- ggplot2::ggplot(
    df_processed,
    ggplot2::aes(
      x = .data[["deg_count"]],
      y = .data[["celltype"]],
      fill = .data[["regulation"]]
    )
  ) +
    ggplot2::geom_bar(
      stat = "identity",
      position = ggplot2::position_stack(
        reverse = TRUE
      )
    ) +
    ggplot2::geom_text(
      ggplot2::aes(
        label = .data[["segment_label"]]
      ),
      position = ggplot2::position_stack(
        reverse = TRUE,
        vjust = 0.5
      ),
      colour = stack_text_color,
      size = text_size,
      fontface = "bold"
    ) +
    ggplot2::geom_text(
      data = total_label_df,
      ggplot2::aes(
        x = .data[["total_deg"]],
        y = .data[["celltype"]],
        label = .data[["label_text"]]
      ),
      inherit.aes = FALSE,
      hjust = -0.25,
      size = text_size,
      fontface = "bold"
    ) +
    ggplot2::scale_fill_manual(
      values = color_vec
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
      x = x_lab,
      y = y_lab,
      fill = "Regulation",
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
        size = base_font_size + 1
      )
    ) +
    ggplot2::guides(
      fill = ggplot2::guide_legend(
        reverse = FALSE
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
