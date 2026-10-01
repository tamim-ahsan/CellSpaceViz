# ============================================================
# Internal helper functions
# ============================================================

# Internal factor-level helper
.get_factor_levels <- function(x) {

  if (is.factor(x)) {

    levels(x)

  } else {

    sort(
      unique(
        stats::na.omit(x)
      )
    )

  }

}


# Internal colour validation helper
.validate_colors <- function(
    levels,
    color_vec = NULL
) {

  if (length(levels) == 0) {
    stop(
      "'levels' must contain at least one value.",
      call. = FALSE
    )
  }

  if (is.null(color_vec)) {

    color_vec <- scales::hue_pal()(length(levels))
    names(color_vec) <- levels

  } else {

    if (is.null(names(color_vec))) {

      if (length(color_vec) < length(levels)) {

        stop(
          "The supplied 'color_vec' contains fewer colors (",
          length(color_vec),
          ") than groups (",
          length(levels),
          ").",
          call. = FALSE
        )

      }

      color_vec <- color_vec[seq_along(levels)]
      names(color_vec) <- levels

    } else {

      missing_levels <- setdiff(
        levels,
        names(color_vec)
      )

      if (length(missing_levels) > 0) {

        stop(
          "Missing colors in 'color_vec' for group(s): ",
          paste(
            missing_levels,
            collapse = ", "
          ),
          call. = FALSE
        )

      }

      color_vec <- color_vec[levels]

    }

  }

  return(color_vec)

}


# Internal required-column validator
.validate_columns <- function(
    df,
    required_cols
) {

  missing_cols <- setdiff(
    required_cols,
    colnames(df)
  )

  if (length(missing_cols) > 0) {

    stop(
      "Missing column(s) in 'df': ",
      paste(
        missing_cols,
        collapse = ", "
      ),
      call. = FALSE
    )

  }

  invisible(TRUE)

}
